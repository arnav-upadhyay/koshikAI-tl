"""
Cross-N analysis plots for the spindler/rplL elongation sweep.

Reads `data/output/spindler_rplL_elongation_log.csv` (per-codon rows for every run)
plus per-N `spindler_rplL_N{N}_first_summary.csv` files (pool snapshots) and produces
five publication-grade comparison figures into `data/output/figures/`.

Canonical-row provenance: see `data/output/runs_summary.md`.

Run from project root:
    python src/plot_cross_N.py

Outputs (all PNG, 150 dpi):
    figures/cross_N_elongation_profile.png       # v_smooth(codon) overlay, one curve per N
    figures/cross_N_sigma_tau_vs_N.png           # Σ τ and mean v vs N (invariance plot)
    figures/cross_N_slowest_codons.png           # τ at top-K slowest positions, per N
    figures/cross_N_pool_stability.png           # T3 / RIB_active / PRT vs N
    figures/cross_N_wall_time.png                # wall time vs N (compute scaling)
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Configuration: canonical row timestamps (see runs_summary.md)
# ---------------------------------------------------------------------------
CANONICAL_TIMESTAMPS: dict[int, str] = {
    10: "2026-05-27T10:12:49.800",
    20: "2026-05-27T15:37:32.860",
    25: "2026-05-22T23:19:21.121",  # tuning-era F scales — flagged in summary
    30: "2026-05-25T15:19:51.616",  # post-warm-start (Run I')
    35: "2026-05-23T16:43:56.419",  # tuning-era K scales — flagged
    40: "2026-05-24T11:35:35.939",  # tuning-era M scales — flagged
    45: "2026-05-25T22:19:22.412",
    50: "2026-05-26T14:32:38.477",
}

# Subset of canonical N values that use *strict* run-I scales
# (V_IN=20, ARS=20, tRNA=1.7, EF=√2, Rib=2). Used for the invariance fit/band.
STRICT_RUN_I = {10, 20, 30, 45, 50}

# Paths (script-location-aware so it works from anywhere)
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
OUTPUT_DIR = PROJECT_ROOT / "data" / "output"
FIGURES_DIR = OUTPUT_DIR / "figures"
ELONGATION_LOG = OUTPUT_DIR / "spindler_rplL_elongation_log.csv"


# ---------------------------------------------------------------------------
# Plot style
# ---------------------------------------------------------------------------
plt.rcParams.update(
    {
        "figure.dpi": 110,
        "savefig.dpi": 150,
        "savefig.bbox": "tight",
        "font.size": 11,
        "axes.titlesize": 12,
        "axes.labelsize": 11,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "legend.frameon": False,
    }
)


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def load_canonical_codon_rows() -> pd.DataFrame:
    """Return the elongation-log rows for the canonical timestamps only.

    Adds a `regime` column: 'strict' if N ∈ STRICT_RUN_I else 'tuning-era'.
    """
    df = pd.read_csv(ELONGATION_LOG)
    df["timestamp"] = df["timestamp"].astype(str)

    keep = []
    for N, ts in CANONICAL_TIMESTAMPS.items():
        mask = (df["timestamp"] == ts) & (df["N"] == N)
        n_match = int(mask.sum())
        if n_match == 0:
            print(f"WARNING: no rows for N={N} ts={ts}")
            continue
        if n_match not in (119, 120, 121):  # rplL has 120 codons; allow ±1 for indexing
            print(f"WARNING: N={N} ts={ts} has {n_match} rows (expected 120)")
        keep.append(df[mask].copy())

    out = pd.concat(keep, ignore_index=True)
    out["regime"] = np.where(out["N"].isin(STRICT_RUN_I), "strict", "tuning-era")
    return out


def load_canonical_summary(codon_df: pd.DataFrame) -> pd.DataFrame:
    """Pool-level metrics for each canonical N.

    The per-N `spindler_rplL_N{N}_first_summary.csv` files are overwritten by every
    fresh attempt for that N (including failed reruns), so they may be stale. We:
      1. always trust the elongation-log row for `T3_mean`, `T3_min`, `m_bg_each`,
         scales, and `wall_time` (those are timestamped per row);
      2. trust the summary CSV for the pool fields NOT present in the elongation log
         (`RIB_active`, `EGT`, `EFG`, `PRT_primary`) **only when** that summary's
         `status == "ok"` AND its `rplL_total_time_s` matches the elongation log
         within 0.5 % (i.e. the summary CSV is from the same run as the canonical
         row, not from a later overwrite by a failed retry).
    """
    log_first = (
        codon_df.groupby("N", as_index=False)
        .agg(log_sigma_tau=("rplL_total_time_s", "first"),
             log_T3_mean=("T3_mean", "first"),
             log_T3_min=("T3_min", "first"))
    )

    pool_fields = ["RIB_active", "EGT", "EFG", "PRT_primary"]
    rows = []
    for N in CANONICAL_TIMESTAMPS:
        path = OUTPUT_DIR / f"spindler_rplL_N{N}_first_summary.csv"
        log_row = log_first[log_first["N"] == N].iloc[0] if (log_first["N"] == N).any() else None
        if log_row is None:
            continue

        rec = {
            "N": N,
            "T3_mean": log_row["log_T3_mean"],
            "T3_min": log_row["log_T3_min"],
        }
        if not path.exists():
            print(f"WARNING: summary CSV missing for N={N} — pool fields will be NaN")
            for f in pool_fields:
                rec[f] = np.nan
            rows.append(rec)
            continue

        s = pd.read_csv(path).iloc[0]
        canonical = (
            s.get("status") == "ok"
            and pd.notna(s.get("rplL_total_time_s"))
            and abs(float(s["rplL_total_time_s"]) - float(log_row["log_sigma_tau"])) < 0.005 * float(log_row["log_sigma_tau"])
        )
        if not canonical:
            print(
                f"NOTE: summary CSV for N={N} does not match the canonical elongation-log row"
                f" (summary status={s.get('status')}, summary sigma_tau={s.get('rplL_total_time_s')}"
                f" vs log sigma_tau={log_row['log_sigma_tau']:.3f}). Pool fields will be NaN."
            )
            for f in pool_fields:
                rec[f] = np.nan
        else:
            for f in pool_fields:
                rec[f] = float(s[f]) if pd.notna(s[f]) else np.nan

        rows.append(rec)

    summary = pd.DataFrame(rows).reset_index(drop=True)
    summary["N"] = summary["N"].astype(int)
    summary["regime"] = np.where(summary["N"].isin(STRICT_RUN_I), "strict", "tuning-era")
    return summary


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def colour_for_N(N_values: list[int]) -> dict[int, tuple]:
    """Map each N to a viridis colour (low N = purple, high N = yellow)."""
    cmap = plt.get_cmap("viridis")
    Ns_sorted = sorted(N_values)
    if len(Ns_sorted) == 1:
        return {Ns_sorted[0]: cmap(0.5)}
    norm = lambda n: (n - Ns_sorted[0]) / (Ns_sorted[-1] - Ns_sorted[0])
    return {n: cmap(norm(n)) for n in Ns_sorted}


def per_N_summary_from_codon_log(codon_df: pd.DataFrame) -> pd.DataFrame:
    """Collapse the codon-level dataframe to one row per N (any column with stable
    value over the 120 codons can be pulled off the first row)."""
    out = (
        codon_df.groupby("N", as_index=False)
        .agg(
            timestamp=("timestamp", "first"),
            sigma_tau=("rplL_total_time_s", "first"),
            mean_v=("rplL_mean_v_cs", "first"),
            T3_mean=("T3_mean", "first"),
            T3_min=("T3_min", "first"),
            wall_time=("wall_time_s", "first"),
            V_IN=("V_IN_70SIC_scale", "first"),
            tRNA=("tRNA_scale", "first"),
            atol=("ode_abstol", "first"),
            m_bg=("m_bg_each", "first"),
        )
        .sort_values("N")
        .reset_index(drop=True)
    )
    out["regime"] = np.where(out["N"].isin(STRICT_RUN_I), "strict", "tuning-era")
    return out


# ---------------------------------------------------------------------------
# Plot 1: smoothed elongation rate per codon — overlay all N
# ---------------------------------------------------------------------------
def plot_elongation_profile(codon_df: pd.DataFrame, out_path: Path) -> None:
    fig, axes = plt.subplots(
        2, 1, figsize=(11, 8), sharex=True, gridspec_kw={"height_ratios": [3, 2]}
    )
    ax_v, ax_tau = axes
    Ns = sorted(codon_df["N"].unique())
    colours = colour_for_N(Ns)

    for N in Ns:
        sub = codon_df[codon_df["N"] == N].sort_values("codon_pos")
        strict = N in STRICT_RUN_I
        ls = "-" if strict else "--"
        lw = 1.5 if strict else 1.2
        alpha = 0.95 if strict else 0.6
        ax_v.plot(
            sub["codon_pos"],
            sub["v_smooth_cod_s"],
            color=colours[N],
            ls=ls,
            lw=lw,
            alpha=alpha,
            label=f"N = {N}" + ("" if strict else "  (tuning-era)"),
        )
        ax_tau.plot(
            sub["codon_pos"],
            sub["tau_s"],
            color=colours[N],
            ls=ls,
            lw=lw,
            alpha=alpha,
        )

    ax_v.set_ylabel("v_smooth  (codons / s)")
    ax_v.set_title(
        "Smoothed rplL elongation-rate profile across background-gene count N\n"
        "(curves should overlap if cross-N invariance holds — sequence-determined, not pool-determined)"
    )
    ax_v.legend(ncol=2, fontsize=9, loc="lower right")
    ax_v.set_ylim(bottom=0)

    ax_tau.set_xlabel("codon position on rplL")
    ax_tau.set_ylabel("τ  (s per codon)")
    ax_tau.set_yscale("log")
    ax_tau.axhline(0.2, color="grey", lw=0.6, ls=":", alpha=0.6)
    ax_tau.text(
        2, 0.21, "τ = 0.2 s  (slow-codon threshold)", color="grey", fontsize=8, va="bottom"
    )

    fig.savefig(out_path)
    plt.close(fig)
    print(f"wrote {out_path}")


# ---------------------------------------------------------------------------
# Plot 2: Σ τ and mean v vs N (invariance metric)
# ---------------------------------------------------------------------------
def plot_sigma_tau_vs_N(per_N: pd.DataFrame, out_path: Path) -> None:
    fig, (ax_s, ax_v) = plt.subplots(1, 2, figsize=(12, 4.5))

    strict_mask = per_N["regime"] == "strict"
    tuning_mask = ~strict_mask

    def scatter(ax, x, y, mask, **kw):
        ax.scatter(x[mask], y[mask], **kw)

    # Σ τ
    scatter(
        ax_s, per_N["N"], per_N["sigma_tau"], strict_mask,
        color="#1f77b4", s=70, zorder=3, label="strict run-I scales",
    )
    scatter(
        ax_s, per_N["N"], per_N["sigma_tau"], tuning_mask,
        facecolor="white", edgecolor="#1f77b4", s=70, zorder=3,
        label="tuning-era scales", linewidth=1.4,
    )

    # band over strict-regime mean ± std
    strict_vals = per_N.loc[strict_mask, "sigma_tau"]
    if len(strict_vals) >= 3:
        m, s = strict_vals.mean(), strict_vals.std(ddof=1)
        ax_s.axhline(m, color="#1f77b4", lw=0.8, alpha=0.6)
        ax_s.fill_between(
            [per_N["N"].min() - 2, per_N["N"].max() + 2],
            m - s, m + s,
            color="#1f77b4", alpha=0.1, zorder=1,
            label=f"strict-regime mean ± std  = {m:.3f} ± {s:.3f} s",
        )

    ax_s.set_xlabel("N  (number of background genes)")
    ax_s.set_ylabel("Σ τ on rplL  (s)")
    ax_s.set_title("rplL total elongation time vs N")
    ax_s.legend(fontsize=9, loc="lower right")
    ax_s.set_ylim(7.30, 7.85)
    ax_s.set_xlim(per_N["N"].min() - 2, per_N["N"].max() + 2)

    # mean v
    scatter(
        ax_v, per_N["N"], per_N["mean_v"], strict_mask,
        color="#d62728", s=70, zorder=3, label="strict run-I",
    )
    scatter(
        ax_v, per_N["N"], per_N["mean_v"], tuning_mask,
        facecolor="white", edgecolor="#d62728", s=70, zorder=3,
        label="tuning-era", linewidth=1.4,
    )
    if len(strict_vals) >= 3:
        mv = per_N.loc[strict_mask, "mean_v"]
        ax_v.axhline(mv.mean(), color="#d62728", lw=0.8, alpha=0.6)
        ax_v.fill_between(
            [per_N["N"].min() - 2, per_N["N"].max() + 2],
            mv.mean() - mv.std(ddof=1), mv.mean() + mv.std(ddof=1),
            color="#d62728", alpha=0.1, zorder=1,
            label=f"mean ± std = {mv.mean():.2f} ± {mv.std(ddof=1):.2f} cod/s",
        )
    ax_v.set_xlabel("N  (number of background genes)")
    ax_v.set_ylabel("rplL mean v  (codons / s)")
    ax_v.set_title("rplL mean elongation rate vs N")
    ax_v.legend(fontsize=9, loc="upper right")
    ax_v.set_ylim(16.5, 18.0)
    ax_v.set_xlim(per_N["N"].min() - 2, per_N["N"].max() + 2)

    fig.suptitle(
        "Cross-N invariance: rplL elongation observables should be ~flat in N",
        y=1.02, fontsize=12,
    )
    fig.savefig(out_path)
    plt.close(fig)
    print(f"wrote {out_path}")


# ---------------------------------------------------------------------------
# Plot 3: limiting / slowest codons — τ at the top-K slowest positions, per N
# ---------------------------------------------------------------------------
def plot_slowest_codons(codon_df: pd.DataFrame, out_path: Path, top_k: int = 12) -> None:
    # Rank slowest positions using the *median* τ across strict-regime runs
    strict_df = codon_df[codon_df["regime"] == "strict"]
    if len(strict_df) == 0:
        strict_df = codon_df  # fallback
    median_tau = (
        strict_df.groupby("codon_pos")["tau_s"]
        .median()
        .sort_values(ascending=False)
    )
    top_positions = median_tau.head(top_k).index.tolist()

    Ns = sorted(codon_df["N"].unique())
    colours = colour_for_N(Ns)

    fig, ax = plt.subplots(figsize=(12, 5.5))
    width = 0.8 / len(Ns)
    x = np.arange(len(top_positions))

    for i, N in enumerate(Ns):
        sub = codon_df[(codon_df["N"] == N) & (codon_df["codon_pos"].isin(top_positions))]
        sub = sub.set_index("codon_pos").reindex(top_positions)["tau_s"]
        offset = (i - (len(Ns) - 1) / 2) * width
        ax.bar(
            x + offset, sub.values, width=width,
            color=colours[N], label=f"N={N}",
            edgecolor="white", linewidth=0.5,
        )

    ax.set_xticks(x)
    ax.set_xticklabels([f"pos {p}" for p in top_positions], rotation=45, ha="right")
    ax.set_ylabel("τ  (s)")
    ax.set_title(
        f"Top {top_k} slowest rplL codons across N\n"
        "(position identity is determined by the rplL sequence; height should be roughly N-invariant)"
    )
    ax.legend(ncol=4, fontsize=9, loc="upper right")
    ax.set_ylim(0, max(median_tau.head(top_k)) * 1.4)

    # annotate which mean-τ ranks these are
    for i, pos in enumerate(top_positions):
        ax.text(
            x[i], median_tau.iloc[i] * 1.05,
            f"#{i+1}", ha="center", va="bottom", fontsize=8, color="#444",
        )

    fig.savefig(out_path)
    plt.close(fig)
    print(f"wrote {out_path}")

    # Also write a small CSV with the identity of the top-K positions
    diag_path = out_path.with_suffix(".csv")
    diag = pd.DataFrame({
        "rank": np.arange(1, len(top_positions) + 1),
        "codon_pos": top_positions,
        "median_tau_strict_s": median_tau.head(top_k).values,
    })
    diag.to_csv(diag_path, index=False)
    print(f"wrote {diag_path}  (slowest-codon position table)")


# ---------------------------------------------------------------------------
# Plot 4: pool stability — T3, RIB_active, EFG/EGT, PRT vs N
# ---------------------------------------------------------------------------
def plot_pool_stability(summary: pd.DataFrame, out_path: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
    strict_mask = summary["regime"] == "strict"
    tuning_mask = ~strict_mask

    def two_class_scatter(ax, x, y, color, ylabel, title):
        ax.scatter(x[strict_mask], y[strict_mask], color=color, s=60, zorder=3, label="strict run-I")
        ax.scatter(
            x[tuning_mask], y[tuning_mask],
            facecolor="white", edgecolor=color, s=60, zorder=3,
            linewidth=1.4, label="tuning-era",
        )
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend(fontsize=8, loc="best")

    two_class_scatter(
        axes[0, 0], summary["N"], summary["T3_mean"], "#2ca02c",
        "T3 mean  (µM)", "Mean charged tRNA·EF-Tu·GTP pool",
    )
    axes[0, 0].axhline(0.4, color="grey", ls=":", lw=0.6)
    axes[0, 0].text(summary["N"].min(), 0.45, "KM_T3 = 0.4 µM", color="grey", fontsize=8)

    two_class_scatter(
        axes[0, 1], summary["N"], summary["RIB_active"], "#9467bd",
        "RIB_active  (µM)", "Active-ribosome pool (sum across all genes)",
    )

    two_class_scatter(
        axes[1, 0], summary["N"], summary["EGT"], "#ff7f0e",
        "EGT  (µM)", "EF-Tu·GTP·tRNA total complex pool",
    )

    two_class_scatter(
        axes[1, 1], summary["N"], summary["PRT_primary"], "#1f77b4",
        "PRT(rplL)  (µM)", "Cumulative rplL protein at t=30 min",
    )

    for ax in axes[1, :]:
        ax.set_xlabel("N  (number of background genes)")

    fig.suptitle("Shared-pool sanity checks — should be roughly invariant in N at strict run-I scales", y=1.0)
    fig.savefig(out_path)
    plt.close(fig)
    print(f"wrote {out_path}")


# ---------------------------------------------------------------------------
# Plot 5: wall time vs N (compute scaling)
# ---------------------------------------------------------------------------
def plot_wall_time(per_N: pd.DataFrame, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    strict_mask = per_N["regime"] == "strict"
    tuning_mask = ~strict_mask

    ax.scatter(per_N.loc[strict_mask, "N"], per_N.loc[strict_mask, "wall_time"],
               color="#1f77b4", s=70, zorder=3, label="strict run-I")
    ax.scatter(per_N.loc[tuning_mask, "N"], per_N.loc[tuning_mask, "wall_time"],
               facecolor="white", edgecolor="#1f77b4", s=70, zorder=3,
               linewidth=1.4, label="tuning-era")

    # Annotate atol for each point
    for _, r in per_N.iterrows():
        ax.annotate(
            f"atol={r['atol']:.0e}",
            (r["N"], r["wall_time"]),
            textcoords="offset points", xytext=(7, -3),
            fontsize=8, color="#555",
        )

    ax.set_xlabel("N  (number of background genes)")
    ax.set_ylabel("wall time  (s)")
    ax.set_title("ODE-solve wall time vs N (atol annotated)")
    ax.set_yscale("log")
    ax.legend(loc="upper left", fontsize=9)

    fig.savefig(out_path)
    plt.close(fig)
    print(f"wrote {out_path}")


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def main() -> None:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Reading {ELONGATION_LOG} …")
    codon_df = load_canonical_codon_rows()
    per_N = per_N_summary_from_codon_log(codon_df)

    print("\nCanonical per-N summary:")
    print(per_N.to_string(index=False))

    print(f"\nReading per-N summary CSVs …")
    summary = load_canonical_summary(codon_df)

    print("\nProducing figures …")
    plot_elongation_profile(codon_df, FIGURES_DIR / "cross_N_elongation_profile.png")
    plot_sigma_tau_vs_N(per_N, FIGURES_DIR / "cross_N_sigma_tau_vs_N.png")
    plot_slowest_codons(codon_df, FIGURES_DIR / "cross_N_slowest_codons.png")
    plot_pool_stability(summary, FIGURES_DIR / "cross_N_pool_stability.png")
    plot_wall_time(per_N, FIGURES_DIR / "cross_N_wall_time.png")

    print("\nDone.")


if __name__ == "__main__":
    main()
