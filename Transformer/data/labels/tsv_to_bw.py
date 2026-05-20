"""
Convert signal files → one BigWig per track (G1, ES, MS, LS).

Rice:
  source: TSV (log1p(soft_clip(CPM))) → expm1 → CPM
  output: bw/rice/

Arabidopsis / Zeamay:
  source: raw-count BED matrix → CPM → replicate mean  (no soft_clip, no log1p)
  output: bw/arabidopsis_cpm/  bw/zeamay_cpm/

Usage:
  python tsv_to_bw.py                    # all three
  python tsv_to_bw.py arabidopsis zeamay
"""
import sys
import numpy as np
import pandas as pd
import pyBigWig
from pathlib import Path

LABELS_DIR = Path(__file__).parent
BED_DIR    = Path('/Users/lzz/Downloads/qq')
TRACKS     = ['G1', 'ES', 'MS', 'LS']

# ── Rice: from TSV ────────────────────────────────────────────────────────────
def tsv_to_bw(species: str, tsv: Path, outdir: Path) -> None:
    print(f'\n[{species}] reading {tsv.name} ...')
    df = pd.read_csv(tsv, sep='\t', dtype={'chrom': str})
    print(f'  {len(df):,} bins')

    # TSV is log1p(soft_clip(CPM)) → expm1 → soft_clipped CPM
    for t in TRACKS:
        df[t] = np.expm1(df[t].values.astype(np.float64)).astype(np.float32)

    _write_bw(df, outdir)


# ── Arabidopsis / Zeamay: from raw BED, no clip ───────────────────────────────
BED_CONFIGS = {
    'arabidopsis': {
        'bed':  BED_DIR / 'arabidopsis_128bp_matrix_sorted.bed',
        'cols': ['G1', 'ES_rep1', 'ES_rep2', 'ES_rep3',
                 'MS_rep1', 'MS_rep2', 'MS_rep3',
                 'LS_rep1', 'LS_rep2', 'LS_rep3'],
        'groups': {
            'G1': ['G1'],
            'ES': ['ES_rep1', 'ES_rep2', 'ES_rep3'],
            'MS': ['MS_rep1', 'MS_rep2', 'MS_rep3'],
            'LS': ['LS_rep1', 'LS_rep2', 'LS_rep3'],
        },
        'outdir': LABELS_DIR / 'bw' / 'arabidopsis_cpm',
    },
    'zeamay': {
        'bed':  BED_DIR / 'zeamay2017_128bp_matrix_sorted.bed',
        'cols': ['G1', 'ES_rep1', 'ES_rep2', 'ES_rep3',
                 'MS_rep1', 'MS_rep2', 'MS_rep3',
                 'LS_rep1', 'LS_rep2'],
        'groups': {
            'G1': ['G1'],
            'ES': ['ES_rep1', 'ES_rep2', 'ES_rep3'],
            'MS': ['MS_rep1', 'MS_rep2', 'MS_rep3'],
            'LS': ['LS_rep1', 'LS_rep2'],
        },
        'outdir': LABELS_DIR / 'bw' / 'zeamay_cpm',
    },
}

def bed_to_bw(species: str, cfg: dict) -> None:
    bed, sample_cols, groups, outdir = cfg['bed'], cfg['cols'], cfg['groups'], cfg['outdir']
    col_names = ['chrom', 'start', 'end'] + sample_cols

    print(f'\n[{species}] pass 1: library sizes from {bed.name} ...')
    lib = np.zeros(len(sample_cols), dtype=np.float64)
    for chunk in pd.read_csv(bed, sep='\t', header=None, names=col_names,
                             dtype={'chrom': str}, chunksize=200_000, low_memory=False):
        lib += chunk[sample_cols].values.astype(np.float64).sum(axis=0)
    lib = lib.clip(min=1)

    print(f'[{species}] pass 2: streaming CPM → replicate mean → bw ...')
    col_idx = {c: i for i, c in enumerate(sample_cols)}
    phase_idxs = {phase: [col_idx[r] for r in reps] for phase, reps in groups.items()}

    # chrom order and sizes (single scan)
    sizes: dict[str, int] = {}
    order: list[str] = []
    for chunk in pd.read_csv(bed, sep='\t', header=None, names=col_names,
                             usecols=['chrom', 'end'], dtype={'chrom': str},
                             chunksize=200_000, low_memory=False):
        for chrom, end in zip(chunk['chrom'], chunk['end']):
            if chrom not in sizes:
                sizes[chrom] = int(end); order.append(chrom)
            elif int(end) > sizes[chrom]:
                sizes[chrom] = int(end)
    header = [(c, sizes[c]) for c in order]

    outdir.mkdir(parents=True, exist_ok=True)
    bw_handles = {}
    for phase in TRACKS:
        bw = pyBigWig.open(str(outdir / f'{phase}.bw'), 'w')
        bw.addHeader(header)
        bw_handles[phase] = bw
        print(f'  writing {outdir}/{phase}.bw ...')

    n = 0
    for chunk in pd.read_csv(bed, sep='\t', header=None, names=col_names,
                             dtype={'chrom': str}, chunksize=200_000, low_memory=False):
        counts = chunk[sample_cols].values.astype(np.float64) / lib * 1e6  # CPM, no clip
        chroms = chunk['chrom'].tolist()
        starts = chunk['start'].tolist()
        ends   = chunk['end'].tolist()
        for phase, idxs in phase_idxs.items():
            vals = counts[:, idxs].mean(axis=1).astype(np.float32).tolist()
            bw_handles[phase].addEntries(chroms, starts, ends=ends, values=vals)
        n += len(chunk)
        if n % 2_000_000 == 0:
            print(f'  {n:,} bins ...', flush=True)

    for bw in bw_handles.values():
        bw.close()
    print(f'  done → {outdir}')


# ── helpers ───────────────────────────────────────────────────────────────────
def _write_bw(df: pd.DataFrame, outdir: Path) -> None:
    chrom_order = list(dict.fromkeys(df['chrom'].tolist()))
    chrom_sizes = df.groupby('chrom')['end'].max()
    header = [(c, int(chrom_sizes[c])) for c in chrom_order]
    outdir.mkdir(parents=True, exist_ok=True)
    for track in TRACKS:
        out = outdir / f'{track}.bw'
        print(f'  writing {out} ...')
        bw = pyBigWig.open(str(out), 'w')
        bw.addHeader(header)
        for chrom in chrom_order:
            grp = df[df['chrom'] == chrom].sort_values('start')
            bw.addEntries([chrom] * len(grp), grp['start'].tolist(),
                          ends=grp['end'].tolist(), values=grp[track].tolist())
        bw.close()
    print(f'  done → {outdir}')


# ── main ──────────────────────────────────────────────────────────────────────
ALL = ['rice', 'arabidopsis', 'zeamay']
targets = sys.argv[1:] if len(sys.argv) > 1 else ALL
for name in targets:
    if name == 'rice':
        tsv_to_bw('rice', LABELS_DIR / 'rice_128bp_rt_signals.tsv', LABELS_DIR / 'bw' / 'rice')
    elif name in BED_CONFIGS:
        bed_to_bw(name, BED_CONFIGS[name])
    else:
        print(f"unknown species '{name}', choices: {ALL}"); sys.exit(1)

print('\nAll done.')
