#!/usr/bin/env python3
"""Robustness tests for CPP-mutant replication-timing polarization."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import pearsonr, spearmanr


MUTANTS = ["tcx2_1", "tcx2_3", "sol1_5", "sol1_8"]
LABELS = {
    "tcx2_1": "oscpp8-1",
    "tcx2_3": "oscpp8-3",
    "sol1_5": "oscpp11-5",
    "sol1_8": "oscpp11-8",
}
COLORS = {
    "tcx2_1": "#C25759",
    "tcx2_3": "#E69191",
    "sol1_5": "#599CB4",
    "sol1_8": "#92B5CA",
}
SEED = 20260826


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--cpp8-bins", required=True)
    parser.add_argument("--cpp11-bins", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--plot-dir", required=True)
    parser.add_argument("--bootstrap", type=int, default=2000)
    return parser.parse_args()


def wrt(es: pd.Series, ms: pd.Series, ls: pd.Series) -> np.ndarray:
    esv, msv, lsv = es.to_numpy(float), ms.to_numpy(float), ls.to_numpy(float)
    denom = esv + msv + lsv
    return np.divide(0.5 * msv + lsv, denom, out=np.full_like(denom, np.nan), where=denom > 0)


def slope_from_sufficient(n: float, sx: float, sy: float, sxx: float, sxy: float) -> float:
    denom = sxx - sx * sx / n
    if n <= 2 or denom <= 0:
        return np.nan
    return (sxy - sx * sy / n) / denom


def block_bootstrap_slope(data: pd.DataFrame, nboot: int, rng: np.random.Generator) -> tuple[float, float, float]:
    d = data[["block", "wt_wrt", "delta"]].dropna()
    d = d.assign(x2=d["wt_wrt"] ** 2, xy=d["wt_wrt"] * d["delta"])
    agg = d.groupby("block", sort=False).agg(
        n=("wt_wrt", "size"), sx=("wt_wrt", "sum"), sy=("delta", "sum"),
        sxx=("x2", "sum"), sxy=("xy", "sum")
    )
    totals = agg.sum()
    point = slope_from_sufficient(*[float(totals[x]) for x in ["n", "sx", "sy", "sxx", "sxy"]])
    arr = agg[["n", "sx", "sy", "sxx", "sxy"]].to_numpy(float)
    nb = len(arr)
    boots = np.empty(nboot)
    for i in range(nboot):
        total = arr[rng.integers(0, nb, nb)].sum(axis=0)
        boots[i] = slope_from_sufficient(*total)
    lo, hi = np.nanpercentile(boots, [2.5, 97.5])
    return point, float(lo), float(hi)


def block_bootstrap_corr(data: pd.DataFrame, xcol: str, ycol: str, nboot: int, rng: np.random.Generator) -> tuple[float, float, float]:
    d = data[["block", xcol, ycol]].dropna().rename(columns={xcol: "x", ycol: "y"})
    d = d.assign(x2=d.x**2, y2=d.y**2, xy=d.x * d.y)
    agg = d.groupby("block", sort=False).agg(
        n=("x", "size"), sx=("x", "sum"), sy=("y", "sum"),
        sxx=("x2", "sum"), syy=("y2", "sum"), sxy=("xy", "sum")
    )

    def corr(v: np.ndarray) -> float:
        n, sx, sy, sxx, syy, sxy = v
        vx = sxx - sx * sx / n
        vy = syy - sy * sy / n
        if n <= 2 or vx <= 0 or vy <= 0:
            return np.nan
        return (sxy - sx * sy / n) / np.sqrt(vx * vy)

    arr = agg[["n", "sx", "sy", "sxx", "syy", "sxy"]].to_numpy(float)
    point = corr(arr.sum(axis=0))
    nb = len(arr)
    boots = np.empty(nboot)
    for i in range(nboot):
        boots[i] = corr(arr[rng.integers(0, nb, nb)].sum(axis=0))
    lo, hi = np.nanpercentile(boots, [2.5, 97.5])
    return float(point), float(lo), float(hi)


def block_bootstrap_quintiles(data: pd.DataFrame, nboot: int, rng: np.random.Generator) -> pd.DataFrame:
    levels = [1, 2, 3, 4, 5]
    blocks = pd.Index(data["block"].drop_duplicates())
    bmap = {b: i for i, b in enumerate(blocks)}
    sums = np.zeros((len(blocks), len(levels)), dtype=float)
    nums = np.zeros_like(sums)
    grouped = data.dropna(subset=["delta", "wt_quintile"]).groupby(["block", "wt_quintile"])["delta"].agg(["sum", "count"])
    for (block, q), row in grouped.iterrows():
        sums[bmap[block], int(q) - 1] = row["sum"]
        nums[bmap[block], int(q) - 1] = row["count"]
    point = sums.sum(axis=0) / nums.sum(axis=0)
    boots = np.empty((nboot, len(levels)))
    nb = len(blocks)
    for i in range(nboot):
        idx = rng.integers(0, nb, nb)
        boots[i] = sums[idx].sum(axis=0) / nums[idx].sum(axis=0)
    lo, hi = np.nanpercentile(boots, [2.5, 97.5], axis=0)
    return pd.DataFrame({"wt_quintile": levels, "mean_delta": point, "ci_low": lo, "ci_high": hi})


def load_bound(path: str, name: str) -> pd.DataFrame:
    x = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2], names=["chr", "start", "end"])
    x = x.drop_duplicates()
    x[name] = 1
    return x


def main() -> None:
    args = parse_args()
    outdir, plotdir = Path(args.output_dir), Path(args.plot_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    plotdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    raw_names = [
        "chr", "start", "end",
        "WT-1-G1", "WT-2-G1", "WT-1-ES", "WT-2-ES", "WT-1-MS", "WT-1-LS",
        "sol1_5-1-G1", "sol1_5-2-G1", "sol1_5-1-ES", "sol1_5-2-ES", "sol1_5-1-MS", "sol1_5-1-LS",
        "sol1_8-1-G1", "sol1_8-2-G1", "sol1_8-1-ES", "sol1_8-2-ES", "sol1_8-1-MS", "sol1_8-1-LS",
        "tcx2_1-1-G1", "tcx2_1-2-G1", "tcx2_1-1-ES", "tcx2_1-2-ES", "tcx2_1-1-MS", "tcx2_1-1-LS",
        "tcx2_3-1-G1", "tcx2_3-1-ES", "tcx2_3-1-MS", "tcx2_3-1-LS",
    ]
    df = pd.read_csv(args.input, sep="\t", header=None, names=raw_names)
    df = df[df["chr"].astype(str).str.fullmatch(r"chr(0[1-9]|1[0-2])")].copy()
    df["start"] = pd.to_numeric(df["start"], errors="coerce").astype("Int64")
    df["end"] = pd.to_numeric(df["end"], errors="coerce").astype("Int64")
    df = df.dropna(subset=["start", "end"])
    df[["start", "end"]] = df[["start", "end"]].astype(int)
    library_cols = raw_names[3:]
    df[library_cols] = df[library_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    length_kb = (df["end"] - df["start"] + 1).to_numpy(float) / 1000.0
    rpk = df[library_cols].to_numpy(float) / length_kb[:, None]
    df[library_cols] = rpk / (rpk.sum(axis=0, keepdims=True) / 1_000_000.0)
    df["block"] = df["chr"].astype(str) + ":" + (df["start"] // 1_000_000).astype(str)

    replicate_map = {
        "WT": {"G1": ["WT-1-G1", "WT-2-G1"], "ES": ["WT-1-ES", "WT-2-ES"], "MS": ["WT-1-MS"], "LS": ["WT-1-LS"]},
        "sol1_5": {"G1": ["sol1_5-1-G1", "sol1_5-2-G1"], "ES": ["sol1_5-1-ES", "sol1_5-2-ES"], "MS": ["sol1_5-1-MS"], "LS": ["sol1_5-1-LS"]},
        "sol1_8": {"G1": ["sol1_8-1-G1", "sol1_8-2-G1"], "ES": ["sol1_8-1-ES", "sol1_8-2-ES"], "MS": ["sol1_8-1-MS"], "LS": ["sol1_8-1-LS"]},
        "tcx2_1": {"G1": ["tcx2_1-1-G1", "tcx2_1-2-G1"], "ES": ["tcx2_1-1-ES", "tcx2_1-2-ES"], "MS": ["tcx2_1-1-MS"], "LS": ["tcx2_1-1-LS"]},
        "tcx2_3": {"G1": ["tcx2_3-1-G1"], "ES": ["tcx2_3-1-ES"], "MS": ["tcx2_3-1-MS"], "LS": ["tcx2_3-1-LS"]},
    }
    for sample, stage_map in replicate_map.items():
        g1 = df[stage_map["G1"]].mean(axis=1)
        stage_norm = {}
        for stage in ["ES", "MS", "LS"]:
            stage_norm[stage] = df[stage_map[stage]].mean(axis=1) / (g1 + 1e-6)
        df[f"{sample}_wrt"] = wrt(stage_norm["ES"], stage_norm["MS"], stage_norm["LS"])
    df["wt_wrt"] = df["WT_wrt"]
    valid_wt = df["wt_wrt"].notna()
    # Rank ties deterministically before qcut to guarantee equal-sized strata.
    ranked_wt = df.loc[valid_wt, "wt_wrt"].rank(method="first")
    df.loc[valid_wt, "wt_quintile"] = pd.qcut(ranked_wt, 5, labels=[1, 2, 3, 4, 5]).astype(int)
    df.loc[valid_wt, "wt_decile"] = pd.qcut(ranked_wt, 10, labels=False).astype(int)
    for mutant in MUTANTS:
        df[f"delta_{mutant}"] = df[f"{mutant}_wrt"] - df["wt_wrt"]

    cpp8 = load_bound(args.cpp8_bins, "bound_cpp8")
    cpp11 = load_bound(args.cpp11_bins, "bound_cpp11")
    df = df.merge(cpp8, on=["chr", "start", "end"], how="left").merge(cpp11, on=["chr", "start", "end"], how="left")
    df[["bound_cpp8", "bound_cpp11"]] = df[["bound_cpp8", "bound_cpp11"]].fillna(0).astype(int)

    effect_rows: list[dict[str, object]] = []
    quintile_frames: list[pd.DataFrame] = []
    for mutant in MUTANTS:
        d = df[["block", "wt_wrt", "wt_quintile", f"delta_{mutant}"]].rename(columns={f"delta_{mutant}": "delta"}).dropna()
        slope, slope_lo, slope_hi = block_bootstrap_slope(d, args.bootstrap, rng)
        model = smf.ols("delta ~ wt_wrt", data=d).fit(cov_type="cluster", cov_kwds={"groups": d["block"]})
        qstats = block_bootstrap_quintiles(d, args.bootstrap, rng)
        qstats["mutant"] = mutant
        qstats["label"] = LABELS[mutant]
        quintile_frames.append(qstats)
        q1, q5 = qstats.iloc[0], qstats.iloc[-1]
        extreme = d[d["wt_quintile"].isin([1, 5])].copy()
        expected = ((extreme["wt_quintile"] == 1) & (extreme["delta"] < 0)) | ((extreme["wt_quintile"] == 5) & (extreme["delta"] > 0))
        effect_rows.append(
            {
                "mutant": mutant,
                "label": LABELS[mutant],
                "n_bins": len(d),
                "n_blocks": d["block"].nunique(),
                "polarization_slope": slope,
                "slope_ci_low": slope_lo,
                "slope_ci_high": slope_hi,
                "cluster_robust_p": float(model.pvalues["wt_wrt"]),
                "q1_mean_delta": q1["mean_delta"],
                "q5_mean_delta": q5["mean_delta"],
                "q5_minus_q1": q5["mean_delta"] - q1["mean_delta"],
                "extreme_direction_concordance": float(expected.mean()),
            }
        )
    effects = pd.DataFrame(effect_rows)
    quintiles = pd.concat(quintile_frames, ignore_index=True)
    effects.to_csv(outdir / "mutant_polarization_effects_003.tsv", sep="\t", index=False)
    quintiles.to_csv(outdir / "mutant_rt_quintile_block_bootstrap_003.tsv", sep="\t", index=False)

    pair_rows: list[dict[str, object]] = []
    for family, a, b in [("OsCPP8", "tcx2_1", "tcx2_3"), ("OsCPP11", "sol1_5", "sol1_8")]:
        xcol, ycol = f"delta_{a}", f"delta_{b}"
        d = df[["block", "wt_quintile", xcol, ycol]].dropna()
        r, lo, hi = block_bootstrap_corr(d, xcol, ycol, args.bootstrap, rng)
        rho, rho_p = spearmanr(d[xcol], d[ycol])
        pair_rows.append(
            {
                "family": family,
                "allele_a": LABELS[a],
                "allele_b": LABELS[b],
                "n_bins": len(d),
                "pearson_r": r,
                "pearson_ci_low": lo,
                "pearson_ci_high": hi,
                "spearman_rho": rho,
                "spearman_p": rho_p,
                "same_sign_fraction": float((np.sign(d[xcol]) == np.sign(d[ycol])).mean()),
            }
        )
    allele_df = pd.DataFrame(pair_rows)
    allele_df.to_csv(outdir / "independent_allele_concordance_003.tsv", sep="\t", index=False)

    bound_rows: list[dict[str, object]] = []
    for factor, bound_col, mutants in [
        ("OsCPP8", "bound_cpp8", ["tcx2_1", "tcx2_3"]),
        ("OsCPP11", "bound_cpp11", ["sol1_5", "sol1_8"]),
    ]:
        for mutant in mutants:
            d = df[["chr", "block", "wt_wrt", "wt_decile", bound_col, f"delta_{mutant}"]].dropna().copy()
            d = d.rename(columns={bound_col: "bound", f"delta_{mutant}": "delta"})
            d["polarization_score"] = d["delta"] * (d["wt_wrt"] - 0.5)
            d["absolute_shift"] = d["delta"].abs()
            for outcome in ["polarization_score", "absolute_shift"]:
                fit = smf.ols(f"{outcome} ~ bound + C(wt_decile) + C(chr)", data=d).fit(
                    cov_type="cluster", cov_kwds={"groups": d["block"]}
                )
                coef = float(fit.params["bound"])
                se = float(fit.bse["bound"])
                bound_rows.append(
                    {
                        "factor": factor,
                        "mutant": LABELS[mutant],
                        "outcome": outcome,
                        "n_bins": len(d),
                        "n_bound": int(d["bound"].sum()),
                        "adjusted_bound_effect": coef,
                        "ci_low": coef - 1.96 * se,
                        "ci_high": coef + 1.96 * se,
                        "cluster_robust_p": float(fit.pvalues["bound"]),
                        "raw_bound_mean": float(d.loc[d.bound == 1, outcome].mean()),
                        "raw_unbound_mean": float(d.loc[d.bound == 0, outcome].mean()),
                    }
                )
    bound_df = pd.DataFrame(bound_rows)
    bound_df.to_csv(outdir / "cuttag_bound_site_adjusted_effects_003.tsv", sep="\t", index=False)

    derived_cols = ["chr", "start", "end", "block", "wt_wrt", "wt_quintile", "bound_cpp8", "bound_cpp11"]
    derived_cols += [f"delta_{m}" for m in MUTANTS]
    df[derived_cols].to_csv(outdir / "mutant_rt_bin_metrics_003.tsv.gz", sep="\t", index=False)

    plt.rcParams.update({"font.size": 8, "pdf.fonttype": 42, "svg.fonttype": "none"})
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))

    for mutant in MUTANTS:
        q = quintiles[quintiles.mutant == mutant]
        axes[0, 0].plot(q.wt_quintile, q.mean_delta, marker="o", color=COLORS[mutant], label=LABELS[mutant])
        axes[0, 0].fill_between(q.wt_quintile, q.ci_low, q.ci_high, color=COLORS[mutant], alpha=0.18)
    axes[0, 0].axhline(0, color="black", lw=0.7, ls="--")
    axes[0, 0].set_xticks([1, 2, 3, 4, 5], ["Q1\nearliest", "Q2", "Q3", "Q4", "Q5\nlatest"])
    axes[0, 0].set_ylabel("ΔWRT (mutant − WT)")
    axes[0, 0].set_title("a  RT redistribution with block-bootstrap CI", loc="left", fontweight="bold")
    axes[0, 0].legend(ncol=2, frameon=False)

    order = list(reversed(MUTANTS))
    eplot = effects.set_index("mutant").loc[order]
    y = np.arange(len(order))
    axes[0, 1].errorbar(
        eplot.polarization_slope, y,
        xerr=[eplot.polarization_slope - eplot.slope_ci_low, eplot.slope_ci_high - eplot.polarization_slope],
        fmt="o", color="#333333", ecolor="#777777", capsize=3,
    )
    axes[0, 1].axvline(0, color="black", lw=0.7, ls="--")
    axes[0, 1].set_yticks(y, [LABELS[x] for x in order])
    axes[0, 1].set_xlabel("Slope: ΔWRT ~ WT WRT")
    axes[0, 1].set_title("b  Positive slope quantifies polarization", loc="left", fontweight="bold")

    sample_rng = np.random.default_rng(SEED)
    for family, a, b, ax_color in [
        ("OsCPP8", "tcx2_1", "tcx2_3", "Reds"),
        ("OsCPP11", "sol1_5", "sol1_8", "Blues"),
    ]:
        d = df[[f"delta_{a}", f"delta_{b}"]].dropna()
        if len(d) > 100_000:
            d = d.iloc[sample_rng.choice(len(d), 100_000, replace=False)]
        axes[1, 0].hexbin(d.iloc[:, 0], d.iloc[:, 1], gridsize=65, mincnt=1, bins="log", cmap=ax_color, alpha=0.65)
    lim = 0.65
    axes[1, 0].plot([-lim, lim], [-lim, lim], color="black", lw=0.7, ls="--")
    axes[1, 0].set_xlim(-lim, lim)
    axes[1, 0].set_ylim(-lim, lim)
    axes[1, 0].set_xlabel("Allele 1/5 ΔWRT")
    axes[1, 0].set_ylabel("Allele 3/8 ΔWRT")
    axes[1, 0].set_title("c  Independent-allele concordance (red CPP8; blue CPP11)", loc="left", fontweight="bold")

    bf = bound_df[bound_df.outcome == "polarization_score"].copy()
    bf["display"] = bf["factor"] + " / " + bf["mutant"]
    bf = bf.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(bf))
    axes[1, 1].errorbar(
        bf.adjusted_bound_effect, y,
        xerr=[bf.adjusted_bound_effect - bf.ci_low, bf.ci_high - bf.adjusted_bound_effect],
        fmt="o", color="#333333", ecolor="#777777", capsize=3,
    )
    axes[1, 1].axvline(0, color="black", lw=0.7, ls="--")
    axes[1, 1].set_yticks(y, bf.display)
    axes[1, 1].set_xlabel("Adjusted bound-site effect on polarization score")
    axes[1, 1].set_title("d  CUT&Tag-bound bins vs matched genomic background", loc="left", fontweight="bold")

    for ax in axes.flat:
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(plotdir / "mutant_rt_robustness_003.pdf", bbox_inches="tight")
    fig.savefig(plotdir / "mutant_rt_robustness_003.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    print("Polarization effects")
    print(effects.to_string(index=False))
    print("\nIndependent allele concordance")
    print(allele_df.to_string(index=False))
    print("\nCUT&Tag-bound adjusted effects")
    print(bound_df.to_string(index=False))


if __name__ == "__main__":
    main()
