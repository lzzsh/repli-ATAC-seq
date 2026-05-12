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


def _make_loaders(species_configs, split: str, cfg: dict,
                  ddp: bool = False, rank: int = 0, world_size: int = 1) -> dict:
    """Create one DataLoader per species for the given split."""
    loaders = {}
    for sp in species_configs:
        ds = RepliSeqDataset(
            [sp], split,
            window_size=cfg["data"]["input_window_length"],
            rc_prob=cfg["augmentation"]["rc_prob"] if split == "train" else 0.0,
            shift_max=cfg["augmentation"].get("shift_max", 0) if split == "train" else 0,
        )
        if ddp:
            sampler = DistributedSampler(ds, num_replicas=world_size, rank=rank,
                                         shuffle=(split == "train"))
            loader = DataLoader(ds, batch_size=cfg["training"]["batch_size"],
                                sampler=sampler, num_workers=4, pin_memory=True)
        else:
            loader = DataLoader(ds, batch_size=cfg["training"]["batch_size"],
                                shuffle=(split == "train"), num_workers=4, pin_memory=True)
        loaders[sp.name] = loader
    return loaders


def _validate(model, val_loaders: dict, criterion, device) -> dict:
    model.eval()
    all_logits, all_labels = [], []
    total_loss = rt_loss = 0.0
    n_batches = 0
    per_species: dict[str, dict] = {}
    with torch.no_grad():
        for sp_name, loader in val_loaders.items():
            sp_logits, sp_labels = [], []
            for batch in loader:
                batch = {k: v.to(device) for k, v in batch.items()}
                raw = model.module if isinstance(model, DDP) else model
                out = raw(batch["one_hot"], head=sp_name)
                losses = criterion(out, batch)
                total_loss += losses["total"].item()
                rt_loss    += losses["rt"].item()
                n_batches  += 1
                sp_logits.append(out["rt_logits"].float().cpu().numpy())
                sp_labels.append(batch["rt_labels"].cpu().numpy())
            sp_log = np.concatenate(sp_logits)
            sp_lab = np.concatenate(sp_labels)
            # exclude ignore bins from per-species metrics
            valid = sp_lab.reshape(-1) != -1
            if valid.any():
                per_species[sp_name] = evaluate_predictions(
                    sp_log.reshape(-1, 4)[valid], sp_lab.reshape(-1)[valid]
                )
            all_logits.append(sp_log)
            all_labels.append(sp_lab)

    # aggregate across species (ignore bins excluded)
    flat_log = np.concatenate(all_logits).reshape(-1, 4)
    flat_lab = np.concatenate(all_labels).reshape(-1)
    valid = flat_lab != -1
    metrics = evaluate_predictions(flat_log[valid], flat_lab[valid]) if valid.any() else {}

    metrics["val_loss_total"] = total_loss / max(n_batches, 1)
    metrics["val_loss_rt"]    = rt_loss    / max(n_batches, 1)

    # per-species macro F1 (used for early stopping)
    if per_species:
        metrics["val_macro_f1_mean"] = float(
            np.mean([m["macro_f1"] for m in per_species.values()])
        )
    for sp_name, m in per_species.items():
        for k, v in m.items():
            if isinstance(v, float):
                metrics[f"{sp_name}/{k}"] = v

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

    world_size = dist.get_world_size() if ddp else 1
    train_loaders = _make_loaders(species_configs, "train", cfg, ddp, rank, world_size)
    val_loaders   = _make_loaders(species_configs, "val",   cfg, ddp, rank, world_size)

    model = Basenji2Model(
        species_configs=species_configs,
        bn_momentum=cfg["model"]["bn_momentum"],
    ).to(device)

    if ddp:
        model = DDP(model, device_ids=[local_rank], find_unused_parameters=True)

    class_weights = cfg.get("loss", {}).get("class_weights", None)
    criterion = RTClassLoss(class_weights=class_weights).to(device)

    optimizer = torch.optim.Adam(
        model.parameters(),
        lr=cfg["training"]["learning_rate"],
        betas=(0.9, 0.999),
    )
    scaler = GradScaler(enabled=cfg["training"]["mixed_precision"])

    start_epoch = 0
    best_score, best_f1, patience, global_step = 0.0, 0.0, 0, 0

    if resume:
        ckpt = torch.load(resume, map_location=device)
        raw = model.module if isinstance(model, DDP) else model
        raw.load_state_dict(ckpt["model"])
        if "optimizer" in ckpt:
            optimizer.load_state_dict(ckpt["optimizer"])
        if "scaler" in ckpt:
            scaler.load_state_dict(ckpt["scaler"])
        start_epoch   = ckpt.get("epoch", 0)
        best_score    = ckpt.get("best_score", 0.0)
        best_f1       = ckpt.get("best_f1", 0.0)
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
    # actual steps per epoch ≈ sum of all per-species loader lengths
    steps_per_epoch = sum(len(l) for l in train_loaders.values())
    total_steps = cfg["training"]["max_epochs"] * steps_per_epoch // accum

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
            for loader in train_loaders.values():
                loader.sampler.set_epoch(epoch)
        model.train()
        optimizer.zero_grad()
        epoch_loss_total = epoch_loss_rt = 0.0
        epoch_batches = 0
        t0 = time.time()

        sp_names = list(train_loaders.keys())
        iterators = {sp: iter(loader) for sp, loader in train_loaders.items()}
        # Each species runs a full epoch; total steps = sum of all loader lengths.
        # Round-robin ensures interleaving: rice, maize, arabidopsis, rice, ...
        # When a species is exhausted it wraps around so others can finish.
        active = set(sp_names)
        step_i = 0
        accum_count = 0
        while active:
            sp_name = sp_names[step_i % len(sp_names)]
            step_i += 1
            if sp_name not in active:
                continue
            try:
                batch = next(iterators[sp_name])
            except StopIteration:
                active.discard(sp_name)
                continue

            batch = {k: v.to(device) for k, v in batch.items()}
            with autocast(enabled=cfg["training"]["mixed_precision"]):
                out = model(batch["one_hot"], head=sp_name)
                losses = criterion(out, batch)
            scaler.scale(losses["total"] / accum).backward()
            epoch_loss_total += losses["total"].item()
            epoch_loss_rt    += losses["rt"].item()
            epoch_batches    += 1
            global_step      += 1
            accum_count      += 1
            if accum_count % accum == 0:
                scaler.unscale_(optimizer)
                torch.nn.utils.clip_grad_norm_(model.parameters(), cfg["training"]["gradient_clip_norm"])
                scaler.step(optimizer)
                scaler.update()
                scheduler.step()
                optimizer.zero_grad()
                accum_count = 0

        # flush remaining gradients at epoch end
        if accum_count > 0:
            scaler.unscale_(optimizer)
            torch.nn.utils.clip_grad_norm_(model.parameters(), cfg["training"]["gradient_clip_norm"])
            scaler.step(optimizer)
            scaler.update()
            scheduler.step()
            optimizer.zero_grad()
            accum_count = 0

        should_stop = torch.tensor(0, device=device)
        if is_master:
            train_loss_total = epoch_loss_total / epoch_batches
            writer.add_scalar("train/loss_total", train_loss_total, epoch + 1)
            writer.add_scalar("train/loss_rt", epoch_loss_rt / epoch_batches, epoch + 1)
            writer.add_scalar("train/lr", optimizer.param_groups[0]["lr"], epoch + 1)

            metrics = _validate(model, val_loaders, criterion, device)
            val_loss = metrics["val_loss_total"]
            val_f1   = metrics.get("val_macro_f1_mean", metrics.get("macro_f1", 0.0))

            writer.add_scalar("val/loss_total", val_loss, epoch + 1)
            writer.add_scalar("val/loss_rt", metrics["val_loss_rt"], epoch + 1)
            for k, v in metrics.items():
                if isinstance(v, float) and not k.startswith("val_loss"):
                    writer.add_scalar(f"val/{k}", v, epoch + 1)
            logger.info(
                f"epoch={epoch+1} time={time.time()-t0:.0f}s "
                f"train_loss={train_loss_total:.4f} val_loss={val_loss:.4f} "
                f"val_macro_f1={val_f1:.4f}"
            )
            if val_f1 > best_score:
                best_score = val_f1
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
                    "best_f1": best_score,
                    "global_step": global_step,
                    "metrics": metrics,
                }, ckpt_dir / "best_model.pt")
                logger.info("  → new best checkpoint saved")
            else:
                patience += 1
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
