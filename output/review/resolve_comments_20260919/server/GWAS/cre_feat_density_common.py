#!/usr/bin/env python
"""Variant density × CRE × feature, COMMON variants only (MAF >= 0.05).

For each (CRE, feature) cell:
  workspace = OCR ∩ feature
  segment = CRE ∩ OCR ∩ feature
  count common variants in segment vs in (workspace - segment)
  Fisher exact test → OR + P
"""
import os, sys, time, re, math
import h5py
import numpy as np
import pandas as pd
from collections import defaultdict
from intervaltree import IntervalTree
from scipy.stats import fisher_exact

H5 = "/public/home/WBXie/RiceVarMap_hzhao/rice4k_all_geno_frequency_res.h5"
OCR_BED = "/public/home/wbxie/work/egwas_cre/filtered/ocr.bed"
DIR = "/public/home/wbxie/work/egwas_cre/filtered"
GFF = "/public/home/jcli/data/reference/represent_gff3/Rice.gff3"
OUT = "/public/home/wbxie/work/egwas_cre/cre_feat_density_common.tsv"
MAF_THRESH = 0.05

# 6 CREs in display order
CRE_BEDS = [
    ("FIMO",         "fimo.bed"),
    ("FT",           "ft.bed"),
    ("HIS-element",  "bsj.bed"),
    ("HIS-positive", "his_positive.bed"),
    ("HIS-negative", "his_negative.bed"),
    ("HIS-variable", "his_variable.bed"),
]
FEATURES = ["Promoter<1kb", "5UTR", "CDS", "3UTR", "Intron", "Promoter1-3kb", "Downstream<1kb", "DistalIntergenic"]

t0 = time.time()

# 1. Load OCR
print(f"[{time.time()-t0:.0f}s] Loading OCR...", flush=True)
ocr_trees = defaultdict(IntervalTree)
ocr_total_bp = 0
with open(OCR_BED) as f:
    for line in f:
        p = line.rstrip().split("\t")
        if len(p) < 3: continue
        s, e = int(p[1]), int(p[2])
        if e > s:
            ocr_trees[p[0]].addi(s, e)
            ocr_total_bp += e - s
print(f"  OCR: {ocr_total_bp:,} bp", flush=True)

# 2. Load CRE trees
print(f"[{time.time()-t0:.0f}s] Loading CRE trees...", flush=True)
cre_trees = {}
for name, fn in CRE_BEDS:
    tr = defaultdict(IntervalTree)
    with open(f"{DIR}/{fn}") as f:
        for line in f:
            p = line.rstrip().split("\t")
            if len(p) < 3: continue
            s, e = int(p[1]), int(p[2])
            if e > s:
                tr[p[0]].addi(s, e)
    cre_trees[name] = tr
    print(f"  {name}: loaded", flush=True)

# 3. Build feature trees from GFF (priority order)
print(f"[{time.time()-t0:.0f}s] Building feature trees from GFF...", flush=True)
utr5 = defaultdict(IntervalTree)
cds = defaultdict(IntervalTree)
utr3 = defaultdict(IntervalTree)
gene_body = defaultdict(IntervalTree)
mrnas = {}

with open(GFF) as f:
    for line in f:
        if line.startswith("#"): continue
        parts = line.rstrip().split("\t")
        if len(parts) < 9: continue
        chrom = parts[0]
        if not re.match(r"^chr(0[1-9]|1[0-2])$", chrom): continue
        ft = parts[2]
        s, e = int(parts[3]), int(parts[4])
        strand = parts[6]
        attrs = parts[8]
        m = re.search(r"ID=([^;]+)", attrs) or re.search(r"Parent=([^;]+)", attrs)
        mrna_id = m.group(1) if m else None
        if ft == "mRNA":
            mrnas[mrna_id] = (chrom, s, e, strand)
        elif ft == "five_prime_UTR":
            utr5[chrom].addi(s-1, e)
        elif ft == "CDS":
            cds[chrom].addi(s-1, e)
        elif ft == "three_prime_UTR":
            utr3[chrom].addi(s-1, e)

prom_short = defaultdict(IntervalTree)
prom_long = defaultdict(IntervalTree)
downstrm = defaultdict(IntervalTree)

for mrna, (chrom, s, e, strand) in mrnas.items():
    gene_body[chrom].addi(s-1, e)
    if strand == "+":
        tss, tts = s, e
        if tss > 1: prom_short[chrom].addi(max(0, tss-1-1000), tss-1)
        if tss > 1001: prom_long[chrom].addi(max(0, tss-1-3000), max(0, tss-1-1000))
        downstrm[chrom].addi(tts, tts+1000)
    else:
        tss, tts = e, s
        prom_short[chrom].addi(tss, tss+1000)
        prom_long[chrom].addi(tss+1000, tss+3000)
        if tts > 1: downstrm[chrom].addi(max(0, tts-1-1000), tts-1)
print(f"  feature trees built, t={time.time()-t0:.0f}s", flush=True)

def classify(chrom, pos):
    """Return feature name (priority order)."""
    s, e = pos-1, pos
    if prom_short[chrom][s:e]: return "Promoter<1kb"
    if utr5[chrom][s:e]: return "5UTR"
    if cds[chrom][s:e]: return "CDS"
    if utr3[chrom][s:e]: return "3UTR"
    if gene_body[chrom][s:e]: return "Intron"
    if prom_long[chrom][s:e]: return "Promoter1-3kb"
    if downstrm[chrom][s:e]: return "Downstream<1kb"
    return "DistalIntergenic"

# 4. Load bp counts for (OCR, feature) and (CRE, feature) from existing variant_density_by_feature.tsv
# That file has all the bp counts we need (computed by previous job)
print(f"[{time.time()-t0:.0f}s] Loading bp counts from variant_density_by_feature.tsv...", flush=True)
bpf = "/public/home/wbxie/work/egwas_cre/variant_density_by_feature.tsv"
bp_map = {}
with open(bpf) as f:
    header = f.readline().rstrip().split("\t")
    col = {h: i for i, h in enumerate(header)}
    for line in f:
        p = line.rstrip().split("\t")
        if len(p) < len(header): continue
        cre = p[col["CRE"]]
        feat = p[col["feature"]]
        cre_feat_bp = int(p[col["CRE_feat_bp"]])
        ocr_feat_bp = int(p[col["OCR_feat_bp"]])
        bp_map[(cre, feat)] = (cre_feat_bp, ocr_feat_bp)
print(f"  loaded {len(bp_map)} bp records", flush=True)

# 5. Stream HDF5 - for each common variant in OCR, classify feature + check CRE membership
print(f"[{time.time()-t0:.0f}s] Streaming HDF5 (common variants only)...", flush=True)
ocr_var_feat = defaultdict(int)  # feat -> common var count in OCR ∩ feat
cre_var_feat = defaultdict(lambda: defaultdict(int))  # cre -> feat -> count

h5 = h5py.File(H5, "r")
all_grp = h5["All"]
chrom_keys = sorted([k for k in all_grp.keys() if k.startswith("chr")])

n_total = 0
n_common = 0
n_in_ocr = 0

for chrom in chrom_keys:
    print(f"  [{time.time()-t0:.0f}s] {chrom}...", flush=True)
    table = all_grp[chrom]["table"]
    n = table.shape[0]
    n_total += n

    CHUNK = 500_000
    for start in range(0, n, CHUNK):
        end = min(start + CHUNK, n)
        sub = table[start:end]
        positions = sub["position"]
        af_major = sub["values_block_2"][:, 0]
        maf = 1.0 - af_major
        maf = np.where(maf > 0.5, 1.0 - maf, maf)
        common_mask = maf >= MAF_THRESH
        n_common += common_mask.sum()
        # Process only common variants
        common_pos = positions[common_mask]
        for pos in common_pos:
            pos = int(pos)
            # OCR check
            if not ocr_trees[chrom][pos:pos+1]:
                continue
            n_in_ocr += 1
            feat = classify(chrom, pos)
            ocr_var_feat[feat] += 1
            for cre_name in cre_trees:
                if cre_trees[cre_name][chrom][pos:pos+1]:
                    cre_var_feat[cre_name][feat] += 1
    print(f"    chr done, t={time.time()-t0:.0f}s", flush=True)

h5.close()
print(f"\n[{time.time()-t0:.0f}s] Done streaming. total={n_total:,} common={n_common:,} in_OCR={n_in_ocr:,}", flush=True)
print(f"OCR common var by feat: {dict(ocr_var_feat)}", flush=True)

# 6. FET per (CRE, feature)
print(f"\n[{time.time()-t0:.0f}s] Computing FET...", flush=True)
print(f"{'CRE':<15}{'feature':<18}{'lead_in_CRE':>12}{'CRE_bp':>12}{'lead_in_OCR':>14}{'OCR_bp':>12}{'OR':>8}{'P':>10}")
print("-"*120)

with open(OUT, "w") as fout:
    fout.write("CRE\tfeature\tn_common_CRE_feat\tCRE_feat_bp\tn_common_OCR_feat\tOCR_feat_bp\tdens_CRE_per_kb\tdens_nonCRE_per_kb\tOR\tCI95_low\tCI95_high\tFET_P\n")
    for cre_name, _ in CRE_BEDS:
        for feat in FEATURES:
            key = (cre_name, feat)
            if key not in bp_map:
                continue
            cre_bp, ocr_bp = bp_map[key]
            a = cre_var_feat[cre_name].get(feat, 0)
            ocr_var = ocr_var_feat.get(feat, 0)
            c = ocr_var - a
            b = cre_bp - a
            d = (ocr_bp - cre_bp) - c
            if min(a, b, c, d) < 1:
                fout.write(f"{cre_name}\t{feat}\t{a}\t{cre_bp}\t{ocr_var}\t{ocr_bp}\tNA\tNA\tNA\tNA\tNA\tNA\n")
                continue
            OR, p = fisher_exact([[a, b], [c, d]], alternative="two-sided")
            se = math.sqrt(1/a + 1/b + 1/c + 1/d)
            log_or = math.log(OR) if OR > 0 else 0
            ci_low = math.exp(log_or - 1.96*se)
            ci_hi  = math.exp(log_or + 1.96*se)
            dens_cre = a / cre_bp * 1000
            dens_non = c / (ocr_bp - cre_bp) * 1000
            print(f"{cre_name:<15}{feat:<18}{a:>12}{cre_bp:>12}{ocr_var:>14}{ocr_bp:>12}{OR:>8.3f}{p:>10.2e}")
            fout.write(f"{cre_name}\t{feat}\t{a}\t{cre_bp}\t{ocr_var}\t{ocr_bp}\t{dens_cre:.4f}\t{dens_non:.4f}\t{OR:.4f}\t{ci_low:.4f}\t{ci_hi:.4f}\t{p:.3e}\n")

print(f"\nSaved {OUT}")
print(f"Total time: {time.time()-t0:.0f}s")
