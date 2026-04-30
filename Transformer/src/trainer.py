import os
import numpy as np
import torch
import torch.nn as nn
import torch.distributed as dist
from torch.nn.parallel import DistributedDataParallel as DDP
from torch.utils.data import DataLoader, DistributedSampler
from torch.cuda.amp import GradScaler, autocast
from torch.utils.tensorboard import SummaryWriter
from pathlib import Path
import yaml
import logging
from tqdm import tqdm

from .data.dataset import RepliSeqDataset, SpeciesBalancedSampler, load_manifest
from .tokenization.tokenizer import KmerTokenizer
from .models.model import DNATransformer, MultiTaskLoss
from .eval import evaluate_predictions, compute_wrt_from_phase_pred

logger = logging.getLogger(__name__)


def _build_class_weights(dataset: RepliSeqDataset) -> torch.Tensor:
    counts = np.zeros(3)
    for s in dataset.samples:
        counts[s["rt_class"]] += 1
    w = 1.0 / np.sqrt(counts + 1e-6)
    return torch.tensor(w / w.sum() * 3, dtype=torch.float32)


def _validate(model, loader, criterion, device) -> dict:
    model.eval()
    pp, pt, wt, cp, ct = [], [], [], [], []
    total_loss = phase_loss = class_loss = 0.0
    n_batches = 0
    with torch.no_grad():
        for batch in loader:
            batch = {k: v.to(device) for k, v in batch.items()}
            out = model.module.forward(batch["input_ids"], batch["species_id"]) \
                if isinstance(model, DDP) else model(batch["input_ids"], batch["species_id"])
            losses = criterion(out, batch)
            total_loss += losses["total"].item()
            phase_loss += losses["phase"].item()
            class_loss += losses["class"].item()
            n_batches += 1
            pp.append(out["phase_pred"].cpu().numpy())
            pt.append(batch["phase_labels"].cpu().numpy())
            wt.append(batch["wrt"].cpu().numpy())
            cp.append(out["class_logits"].argmax(-1).cpu().numpy())
            ct.append(batch["rt_class"].cpu().numpy())
    metrics = evaluate_predictions(
        np.concatenate(pp), np.concatenate(pt),
        np.concatenate(wt), np.concatenate(cp), np.concatenate(ct),
    )
    metrics["val_loss_total"] = total_loss / n_batches
    metrics["val_loss_phase"] = phase_loss / n_batches
    metrics["val_loss_class"] = class_loss / n_batches
    return metrics


def train(config_path: str):
    with open(config_path) as f:
        cfg = yaml.safe_load(f)

    # DDP setup
    ddp = int(os.environ.get("RANK", -1)) != -1
    if ddp:
        dist.init_process_group("nccl")
        rank = dist.get_rank()
        local_rank = int(os.environ["LOCAL_RANK"])
        device = torch.device(f"cuda:{local_rank}")
        torch.cuda.set_device(device)
    else:
        rank = 0
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    is_master = rank == 0
    torch.manual_seed(cfg["training"]["seed"] + rank)

    if is_master:
        logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")

    species_configs = load_manifest(cfg["data"]["manifest"])
    tokenizer = KmerTokenizer(
        k=cfg["tokenizer"]["k"],
        stride=cfg["tokenizer"]["stride"],
        add_cls=cfg["tokenizer"].get("add_cls_token", True),
    )

    train_ds = RepliSeqDataset(
        species_configs, "train", tokenizer,
        window_size=cfg["data"]["input_window_length"],
        rc_prob=cfg["augmentation"]["rc_prob"],
    )
    val_ds = RepliSeqDataset(
        species_configs, "val", tokenizer,
        window_size=cfg["data"]["input_window_length"],
        rc_prob=0.0,
    )

    if ddp:
        train_sampler = DistributedSampler(train_ds, shuffle=True)
        train_loader = DataLoader(train_ds, batch_size=cfg["training"]["batch_size"],
                                  sampler=train_sampler, num_workers=4, pin_memory=True)
    else:
        bal_sampler = SpeciesBalancedSampler(train_ds, cfg["training"]["batch_size"])
        train_loader = DataLoader(train_ds, batch_sampler=bal_sampler, num_workers=4, pin_memory=True)

    val_loader = DataLoader(val_ds, batch_size=cfg["training"]["batch_size"],
                            shuffle=False, num_workers=4, pin_memory=True)

    class_weights = _build_class_weights(train_ds).to(device)
    model = DNATransformer(
        vocab_size=tokenizer.vocab_size,
        n_species=len(species_configs),
        d_model=cfg["model"]["d_model"],
        n_layers=cfg["model"]["n_layers"],
        n_heads=cfg["model"]["n_heads"],
        dim_feedforward=cfg["model"]["dim_feedforward"],
        dropout=cfg["model"]["dropout"],
        attn_dropout=cfg["model"]["attention_dropout"],
        species_emb_dim=cfg["model"]["species_embedding_dim"],
        head_hidden_dim=cfg["model"]["head_hidden_dim"],
        head_dropout=cfg["model"]["head_dropout"],
    ).to(device)

    if ddp:
        model = DDP(model, device_ids=[local_rank])

    criterion = MultiTaskLoss(
        class_weights=class_weights,
        lambda_phase=cfg["loss"]["lambda_phase"],
        lambda_class=cfg["loss"]["lambda_class"],
    )
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=cfg["training"]["learning_rate"],
        weight_decay=cfg["training"]["weight_decay"],
    )
    total_steps = cfg["training"]["max_epochs"] * len(train_loader)
    warmup = cfg["training"]["warmup_steps"]
    scheduler = torch.optim.lr_scheduler.OneCycleLR(
        optimizer, max_lr=cfg["training"]["learning_rate"],
        total_steps=total_steps, pct_start=warmup / total_steps,
        anneal_strategy="cos",
        final_div_factor=cfg["training"]["learning_rate"] / cfg["training"]["min_lr"],
    )
    scaler = GradScaler(enabled=cfg["training"]["mixed_precision"])

    out_dir = Path(cfg.get("output_dir", "outputs"))
    ckpt_dir = out_dir / "checkpoints"
    log_dir = Path(cfg.get("log_dir", "logs")) / Path(config_path).stem
    if is_master:
        ckpt_dir.mkdir(parents=True, exist_ok=True)
        log_dir.mkdir(parents=True, exist_ok=True)
        writer = SummaryWriter(log_dir=str(log_dir))

    best_score, patience, global_step = -float("inf"), 0, 0
    accum = cfg["training"]["gradient_accumulation_steps"]

    for epoch in range(cfg["training"]["max_epochs"]):
        if ddp:
            train_sampler.set_epoch(epoch)
        model.train()
        optimizer.zero_grad()
        epoch_loss_total = epoch_loss_phase = epoch_loss_class = 0.0
        epoch_batches = 0
        for batch in (tqdm(train_loader, desc=f"Epoch {epoch+1}") if is_master else train_loader):
            batch = {k: v.to(device) for k, v in batch.items()}
            with autocast(enabled=cfg["training"]["mixed_precision"]):
                out = model(batch["input_ids"], batch["species_id"])
                losses = criterion(out, batch)
            scaler.scale(losses["total"] / accum).backward()
            epoch_loss_total += losses["total"].item()
            epoch_loss_phase += losses["phase"].item()
            epoch_loss_class += losses["class"].item()
            epoch_batches += 1
            global_step += 1
            if global_step % accum == 0:
                scaler.unscale_(optimizer)
                torch.nn.utils.clip_grad_norm_(model.parameters(), cfg["training"]["gradient_clip_norm"])
                scaler.step(optimizer)
                scaler.update()
                optimizer.zero_grad()
                scheduler.step()

        # end of epoch validation
        should_stop = torch.tensor(0, device=device)
        if is_master:
            train_loss_total = epoch_loss_total / epoch_batches
            train_loss_phase = epoch_loss_phase / epoch_batches
            train_loss_class = epoch_loss_class / epoch_batches
            writer.add_scalar("train/loss_total", train_loss_total, epoch + 1)
            writer.add_scalar("train/loss_phase", train_loss_phase, epoch + 1)
            writer.add_scalar("train/loss_class", train_loss_class, epoch + 1)
            writer.add_scalar("train/lr", scheduler.get_last_lr()[0], epoch + 1)

            metrics = _validate(model, val_loader, criterion, device)
            score = -metrics["val_loss_total"]
            writer.add_scalar("val/loss_total", metrics["val_loss_total"], epoch + 1)
            writer.add_scalar("val/loss_phase", metrics["val_loss_phase"], epoch + 1)
            writer.add_scalar("val/loss_class", metrics["val_loss_class"], epoch + 1)
            for k, v in metrics.items():
                if isinstance(v, float) and not k.startswith("val_loss"):
                    writer.add_scalar(f"val/{k}", v, epoch + 1)
            logger.info(
                f"epoch={epoch+1} "
                f"train_loss={train_loss_total:.4f} "
                f"val_loss={metrics['val_loss_total']:.4f} "
                f"val_phase={metrics['val_loss_phase']:.4f} "
                f"val_class={metrics['val_loss_class']:.4f}"
            )
            if score > best_score:
                best_score = score
                patience = 0
                raw = model.module if isinstance(model, DDP) else model
                torch.save({"model": raw.state_dict(), "cfg": cfg, "epoch": epoch + 1,
                            "metrics": metrics}, ckpt_dir / "best_model.pt")
                logger.info("  → new best checkpoint saved")
            else:
                patience += 1
            if patience >= cfg["training"]["early_stopping_patience"]:
                logger.info("Early stopping.")
                should_stop = torch.tensor(1, device=device)
        if ddp:
            dist.broadcast(should_stop, src=0)
        if should_stop.item():
            if ddp:
                dist.destroy_process_group()
            return
        model.train()

    if ddp:
        dist.destroy_process_group()
    if is_master:
        writer.close()
