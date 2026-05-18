"""
Convert rice_128bp_rt_signals.tsv → one BigWig per track (G1, ES, MS, LS).
Output: data/labels/bw/rice_{track}.bw
"""
import pandas as pd
import pyBigWig
from pathlib import Path

TSV   = Path(__file__).parent / 'rice_128bp_rt_signals.tsv'
OUTDIR = Path(__file__).parent / 'bw'
OUTDIR.mkdir(exist_ok=True)

TRACKS = ['G1', 'ES', 'MS', 'LS']

print('Reading TSV...')
df = pd.read_csv(TSV, sep='\t', dtype={'chrom': str, 'start': int, 'end': int,
                                        'G1': float, 'ES': float, 'MS': float, 'LS': float})

# chrom sizes from data
chrom_sizes = df.groupby('chrom')['end'].max().sort_index()
header = [(c, int(s)) for c, s in chrom_sizes.items()]

for track in TRACKS:
    out = OUTDIR / f'rice_{track}.bw'
    print(f'Writing {out} ...')
    bw = pyBigWig.open(str(out), 'w')
    bw.addHeader(header)

    for chrom, grp in df.groupby('chrom', sort=True):
        grp = grp.sort_values('start')
        chroms = [chrom] * len(grp)
        starts = grp['start'].tolist()
        ends   = grp['end'].tolist()
        vals   = grp[track].tolist()
        bw.addEntries(chroms, starts, ends=ends, values=vals)

    bw.close()
    print(f'  done: {out}')

print('All done.')
