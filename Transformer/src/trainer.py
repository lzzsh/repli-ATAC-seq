import math
import os
import sys
import time
import numpy as np
import torch
import torch.distributed as dist
from torch.nn.parallel import DistributedDataParallel as DDP
from torch.utils.data import DataLoader, DistributedSampler
from torch.cuda.amp import GradScaler, autocast
from torch.utils.tensorboard import SummaryWriter
from pathlib import Path
import yaml
import logging

from .data.dataset import RepliSeqDataset, load_manifest
from .models.model import Basenji2Model, RTClassLoss
from .eval import evaluate_predictions

logger = logging.getLogger(__name__)


def _forward_multi_species(
    model, batch: dict, criterion, id_to_name: dict
) -> tuple[dict, torch.Tensor]:
    """Group batch by species_id, forward each group through its head.
    Returns (weighted-mean losses dict, logits [B, 896, 4] in original order)."""
    one_hot     = batch["one_hot"]
    rt_labels   = batch["rt_labels"]
    species_ids = batch["species_id"]
    B = one_hot.shape[0]

    logits_out = torch.empty(B, 896, 4, device=one_hot.device)
    total_loss = rt_loss = 0.0

    for sp_id in species_ids.unique():
        sp_name = id_to_name[sp_id.item()]
        idx = (species_ids == sp_id).nonzero(as_tuple=True)[0]
        out = model(one_hot[idx], head=sp_name)
        logits_out[idx] = out["rt_logits"]
        sub_batch = {"rt_labels": rt_labels[idx]}
        losses = criterion(out, sub_batch)
        n = idx.shape[0]
        total_loss += losses["total"] * n
        rt_loss    += losses["rt"]    * n

    return {"total": total_loss / B, "rt": rt_loss / B}, logits_out


def _validate(model, loader, criterion, device, id_to_name: dict) -> dict:
    model.eval()
    all_logits, all_labels = [], []
    total_loss = rt_loss = 0.0
    n_batches = 0
    with torch.no_grad():
        for batch in loader:
            one_hot    = batch["one_hot"].to(device)
            rt_labels  = batch["rt_labels"].to(device)
            species_id = batch["species_id"].to(device)
            dev_batch  = {**batch, "one_hot": one_hot, "rt_labels": rt_labels,
                          "species_id": species_id}
            raw = model.module if isinstance(model, DDP) else model
            losses, logits = _forward_multi_species(raw, dev_batch, criterion, id_to_name)
            total_loss += losses["total"].item()
            rt_loss    += losses["rt"].item()
            n_batches  += 1
            all_logits.append(logits.cpu().numpy())
            all_labels.append(rt_labels.cpu().numpy())

    metrics = evaluate_predictions(np.concatenate(all_logits), np.concatenate(all_labels))
    metrics["val_loss_total"] = total_loss / n_batches
    metrics["val_loss_rt"]    = rt_loss    / n_batches
    return metrics


def train(config_path: str, resume: str | None = None):
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
    id_to_name = {sp.species_id: sp.name for sp in species_configs}

    train_ds = RepliSeqDataset(
        species_configs, "train",
        window_size=cfg["data"]["input_window_length"],
        rc_prob=cfg["augmentation"]["rc_prob"],
        shift_max=cfg["augmentation"].get("shift_max", 0),
    )
    val_ds = RepliSeqDataset(
        species_configs, "val",
        window_size=cfg["data"]["input_window_length"],
        rc_prob=0.0,
    )

    if ddp:
        train_sampler = DistributedSampler(train_ds, shuffle=True)
        train_loader = DataLoader(train_ds, batch_size=cfg["training"]["batch_size"],
                                  sampler=train_sampler, num_workers=4, pin_memory=True)
    else:
        train_loader = DataLoader(train_ds, batch_size=cfg["training"]["batch_size"],
                                  shuffle=True, num_workers=4, pin_memory=True)

    val_loader = DataLoader(val_ds, batch_size=cfg["training"]["batch_size"],
                            shuffle=False, num_workers=4, pin_memory=True)

    model = Basenji2Model(
        species_configs=species_configs,
        bn_momentum=cfg["model"]["bn_momentum"],
    ).to(device)

    if ddp:
        model = DDP(model, device_ids=[local_rank])

    class_weights = cfg.get("loss", {}).get("class_weights", None)
    criterion = RTClassLoss(class_weights=class_weights).to(device)

    optimizer = torch.optim.Adam(
        model.parameters(),
        lr=cfg["training"]["learning_rate"],
        betas=(0.9, 0.999),
    )
    scaler = GradScaler(enabled=cfg["training"]["mixed_precision"])

    start_epoch = 0
    best_score, best_f1, patience, global_step = float("inf"), 0.0, 0, 0

    if resume:
        ckpt = torch.load(resume, map_location=device)
        raw = model.module if isinstance(model, DDP) else model
        raw.load_state_dict(ckpt["model"])
        if "optimizer" in ckpt:
            optimizer.load_state_dict(ckpt["optimizer"])
        if "scaler" in ckpt:
            scaler.load_state_dict(ckpt["scaler"])
        start_epoch   = ckpt.get("epoch", 0)
        best_score    = ckpt.get("best_score", float("inf"))
        best_f1       = ckpt.get("best_score", 0.0)
        global_step   = ckpt.get("global_step", 0)
        if is_master:
            logger.info(f"Resumed from {resume} (epoch {start_epoch}, best_score={best_score:.4f})")
    out_dir = Path(cfg.get("output_dir", "outputs"))
    ckpt_dir = out_dir / "checkpoints"
    log_dir = Path(cfg.get("log_dir", "logs")) / Path(config_path).stem
    if is_master:
        ckpt_dir.mkdir(parents=True, exist_ok=True)
        log_dir.mkdir(parents=True, exist_ok=True)
        writer = SummaryWriter(log_dir=str(log_dir))

    accum = cfg["training"]["gradient_accumulation_steps"]
    early_stop_patience = cfg["training"]["early_stopping_patience"]

    warmup_steps = cfg["training"].get("warmup_steps", 1000)
    total_steps = cfg["training"]["max_epochs"] * len(train_loader) // accum

    def _lr_lambda(step: int) -> float:
        if step < warmup_steps:
            return step / max(1, warmup_steps)
        progress = (step - warmup_steps) / max(1, total_steps - warmup_steps)
        return max(0.0, 0.5 * (1.0 + math.cos(math.pi * progress)))

    scheduler = torch.optim.lr_scheduler.LambdaLR(optimizer, _lr_lambda)

    if resume:
        if "scheduler" in ckpt:
            scheduler.load_state_dict(ckpt["scheduler"])

    for epoch in range(start_epoch, cfg["training"]["max_epochs"]):
        if ddp:
            train_sampler.set_epoch(epoch)
        model.train()
        optimizer.zero_grad()
        epoch_loss_total = epoch_loss_rt = 0.0
        epoch_batches = 0
        t0 = time.time()
        for batch in train_loader:
            batch = {k: v.to(device) for k, v in batch.items()}
            with autocast(enabled=cfg["training"]["mixed_precision"]):
                losses, _ = _forward_multi_species(
                    model.module if isinstance(model, DDP) else model,
                    batch, criterion, id_to_name
                )
            scaler.scale(losses["total"] / accum).backward()
            epoch_loss_total += losses["total"].item()
            epoch_loss_rt += losses["rt"].item()
            epoch_batches += 1
            global_step += 1
            if global_step % accum == 0:
                scaler.unscale_(optimizer)
                torch.nn.utils.clip_grad_norm_(model.parameters(), cfg["training"]["gradient_clip_norm"])
                scaler.step(optimizer)
                scaler.update()
                scheduler.step()
                optimizer.zero_grad()

        should_stop = torch.tensor(0, device=device)
        if is_master:
            train_loss_total = epoch_loss_total / epoch_batches
            writer.add_scalar("train/loss_total", train_loss_total, epoch + 1)
            writer.add_scalar("train/loss_rt", epoch_loss_rt / epoch_batches, epoch + 1)
            writer.add_scalar("train/lr", optimizer.param_groups[0]["lr"], epoch + 1)

            metrics = _validate(model, val_loader, criterion, device, id_to_name)
            val_loss = metrics["val_loss_total"]

            writer.add_scalar("val/loss_total", val_loss, epoch + 1)
            writer.add_scalar("val/loss_rt", metrics["val_loss_rt"], epoch + 1)
            for k, v in metrics.items():
                if isinstance(v, float) and not k.startswith("val_loss"):
                    writer.add_scalar(f"val/{k}", v, epoch + 1)
            logger.info(
                f"epoch={epoch+1} time={time.time()-t0:.0f}s "
                f"train_loss={train_loss_total:.4f} val_loss={val_loss:.4f} "
                f"macro_f1={metrics.get('macro_f1', float('nan')):.4f} "
                f"acc_ES={metrics.get('acc_ES', float('nan')):.4f} "
                f"acc_MS={metrics.get('acc_MS', float('nan')):.4f} "
                f"acc_LS={metrics.get('acc_LS', float('nan')):.4f} "
                f"acc_NR={metrics.get('acc_NR', float('nan')):.4f}"
            )
            if val_loss < best_score:
                best_score = val_loss
                patience = 0
                raw = model.module if isinstance(model, DDP) else model
                torch.save({
                    "model": raw.state_dict(),
                    "optimizer": optimizer.state_dict(),
                    "scheduler": scheduler.state_dict(),
                    "scaler": scaler.state_dict(),
                    "cfg": cfg,
                    "epoch": epoch + 1,
                    "best_score": best_score,
                    "global_step": global_step,
                    "metrics": metrics,
                }, ckpt_dir / "best_model.pt")
                logger.info("  → new best checkpoint saved")
            else:
                patience += 1
            macro_f1 = metrics.get("macro_f1", 0.0)
            if macro_f1 > best_f1:
                best_f1 = macro_f1
            if patience >= early_stop_patience:
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
