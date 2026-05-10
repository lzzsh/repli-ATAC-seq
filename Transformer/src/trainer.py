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


def _validate(model, loader, criterion, device) -> dict:
    model.eval()
    pp, pt = [], []
    total_loss = rt_loss = 0.0
    n_batches = 0
    with torch.no_grad():
        for batch in loader:
            batch = {k: v.to(device) for k, v in batch.items()}
            out = model.module.forward(batch["one_hot"]) \
                if isinstance(model, DDP) else model(batch["one_hot"])
            losses = criterion(out, batch)
            total_loss += losses["total"].item()
            rt_loss += losses["rt"].item()
            n_batches += 1
            pp.append(out["rt_logits"].cpu().numpy())
            pt.append(batch["rt_labels"].cpu().numpy())

    metrics = evaluate_predictions(np.concatenate(pp), np.concatenate(pt))
    metrics["val_loss_total"] = total_loss / n_batches
    metrics["val_loss_rt"] = rt_loss / n_batches
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
        bn_momentum=cfg["model"]["bn_momentum"],
    ).to(device)

    if ddp:
        model = DDP(model, device_ids=[local_rank])

    criterion = RTClassLoss()

    optimizer = torch.optim.SGD(
        model.parameters(),
        lr=cfg["training"]["learning_rate"],
        momentum=cfg["training"].get("momentum", 0.99),
        weight_decay=cfg["training"].get("weight_decay", 1e-4),
    )
    # ReduceLROnPlateau: halve lr when val loss stops improving
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer,
        mode="min",
        factor=cfg["training"]["lr_decay_factor"],
        patience=cfg["training"]["lr_patience"],
        min_lr=cfg["training"]["min_lr"],
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
        if "scheduler" in ckpt:
            scheduler.load_state_dict(ckpt["scheduler"])
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
                out = model(batch["one_hot"])
                losses = criterion(out, batch)
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
                optimizer.zero_grad()

        should_stop = torch.tensor(0, device=device)
        if is_master:
            train_loss_total = epoch_loss_total / epoch_batches
            writer.add_scalar("train/loss_total", train_loss_total, epoch + 1)
            writer.add_scalar("train/loss_rt", epoch_loss_rt / epoch_batches, epoch + 1)
            writer.add_scalar("train/lr", optimizer.param_groups[0]["lr"], epoch + 1)

            metrics = _validate(model, val_loader, criterion, device)
            val_loss = metrics["val_loss_total"]
            scheduler.step(val_loss)

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
