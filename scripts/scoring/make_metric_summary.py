#!/usr/bin/env python3
"""
Summarize four key fine-screening metrics from compiled_scores.csv and plot their distributions.

Metrics (lower is better for all four):
  vina_affinity           AutoDock Vina docking score (kcal/mol)
  rfaa_pae_inter          RoseTTAFold All-Atom protein-ligand PAE (A)
  boltz2_affinity         Boltz-2 predicted affinity (log10 IC50, uM)
  af3_chain_pair_pae_min  AlphaFold 3 min protein-ligand PAE (A); -1 sentinel from af3_scores.py = missing
Metrics whose method was not run (column absent) are skipped with a warning.

Usage:
  python make_metric_summary.py --target-dir <TASK_ROOT>/<PROTEIN> --target-name <PROTEIN> \
      [--controls controls.csv] [--scores path/to/compiled_scores.csv] [--output-dir DIR]

Inputs:
  --scores    default <target-dir>/fine_screening/compiled_scores.csv (output of compile_scores.py)
  --controls  optional table with SMILES + the same metric columns; optional `name` column is used
              in the legend (otherwise "control N"). Controls are numbered in file order.

Outputs (default <target-dir>/fine_screening/):
  metric_summary.csv           SMILES + metrics + <metric>_pct per compound
                               (pct = % of screen the compound beats; higher = better)
  metric_summary_controls.csv  same columns for the controls (only with --controls)
  <target-name>_metric_kde.png histogram + KDE per metric, numbered grey lines per control
                               (skipped with a warning if seaborn/matplotlib are not installed)
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# (column, panel title, x-axis label)
METRICS = [
    ("vina_affinity", "AutoDock Vina", "affinity (kcal/mol)"),
    ("rfaa_pae_inter", "RoseTTAFold All-Atom", "PAE inter (Å)"),
    ("boltz2_affinity", "Boltz-2", "affinity (log10 IC50, µM)"),
    ("af3_chain_pair_pae_min", "AlphaFold 3", "chain-pair PAE min (Å)"),
]
SENTINEL = {"af3_chain_pair_pae_min": -1.0}  # af3_scores.py writes -1 when the JSON value is None
FILL, LINE, INK, CTRL = "#c9d9ea", "#2b5d8a", "#222222", "#8a8a8a"
LINE_TOP = 0.76  # control lines run from the x-axis to this fraction of the panel height


def load_metrics(path):
    """Read a scores table; keep SMILES/name and the metric columns present (sentinels -> NaN)."""
    df = pd.read_csv(path, low_memory=False)
    if "SMILES" not in df.columns:
        sys.exit(f"ERROR: {path} has no SMILES column")
    keep = [c for c in ("name", "SMILES") if c in df.columns] + [m for m, _, _ in METRICS if m in df.columns]
    df = df[keep].copy()
    for m, _, _ in METRICS:
        if m in df.columns:
            df[m] = pd.to_numeric(df[m], errors="coerce")
            if m in SENTINEL:
                df[m] = df[m].mask(df[m] == SENTINEL[m])
    return df


def percentile_vs(screen, values):
    """% of the screen's non-missing values that each value beats (all metrics: lower is better)."""
    ref = np.sort(screen.dropna().to_numpy())
    pct = 100.0 * (len(ref) - np.searchsorted(ref, values, side="right")) / len(ref)
    return np.where(np.isnan(values), np.nan, np.round(pct, 1))


def add_percentiles(df, screen, metrics):
    for m, _, _ in metrics:
        df[f"{m}_pct"] = percentile_vs(screen[m], df[m].to_numpy(dtype=float))
    return df


def control_label(row, i):
    return row["name"] if "name" in row and pd.notna(row["name"]) else f"control {i + 1}"


def label_tiers(xs, span, min_gap=0.035, n_tiers=3):
    """Stagger number labels vertically when control lines are closer than min_gap * x-range."""
    tiers, placed = {}, []  # placed: (x, tier)
    for i, x in sorted(xs.items(), key=lambda kv: kv[1]):
        used = {t for px, t in placed if abs(x - px) < min_gap * span}
        tiers[i] = next((t for t in range(n_tiers) if t not in used), len(used) % n_tiers)
        placed.append((x, tiers[i]))
    return tiers


def plot(screen, controls, metrics, out_path):
    """Histogram + KDE per metric; controls as numbered grey dashed lines + side legend."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import seaborn as sns
    except ImportError as exc:
        print(f"WARNING: skipping plot ({exc}); run with an env that has seaborn, e.g. conplex_test")
        return False
    sns.set_theme(style="ticks", context="notebook",
                  rc={"axes.edgecolor": "#444444", "axes.labelcolor": INK, "xtick.color": INK,
                      "ytick.color": INK, "axes.titlesize": 11, "axes.labelsize": 9.5,
                      "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
                      # Arial on macOS; Liberation Sans (Arial-metric clone) on the cluster
                      "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"]})
    has_ctrl = controls is not None and len(controls) > 0
    ncols = 1 if len(metrics) == 1 else 2
    nrows = -(-len(metrics) // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(11 if ncols == 2 else 6.5, 3.5 * nrows), squeeze=False)
    for ax in axes.flat[len(metrics):]:
        ax.axis("off")
    for ax, (col, title, xlabel) in zip(axes.flat, metrics):
        vals = screen[col].dropna()
        # density-scaled histogram (white bar edges for definition) under the KDE line
        sns.histplot(x=vals, ax=ax, stat="density", bins=45, color=FILL, alpha=0.9,
                     edgecolor="white", linewidth=0.5)
        sns.kdeplot(x=vals, ax=ax, color=LINE, linewidth=1.4)
        ax.set_ylim(0, ax.get_ylim()[1] * 1.35)  # headroom: peak sits ~0.74 of the axes height
        if has_ctrl:
            xs = {i: row[col] for i, (_, row) in enumerate(controls.iterrows()) if pd.notna(row[col])}
            lo, hi = ax.get_xlim()
            tiers = label_tiers(xs, hi - lo, min_gap=0.03, n_tiers=4)
            for i, x in xs.items():
                # short line ending just above the data; number stacked above it
                ax.axvline(x, ymin=0, ymax=LINE_TOP, color=CTRL, linestyle="--", linewidth=1.0, zorder=3)
                ax.text(x, LINE_TOP + 0.015 + 0.055 * tiers[i], str(i + 1), transform=ax.get_xaxis_transform(),
                        ha="center", va="bottom", fontsize=8.5, color=INK, zorder=4)
        ax.set_title(f"{title}  (n = {len(vals):,})", loc="left")
        ax.set_xlabel(f"{xlabel}  ← better")
        ax.set_ylabel("")
        ax.set_yticks([])
        sns.despine(ax=ax, left=True)
    if has_ctrl:
        lines = [f"{i + 1}   {control_label(row, i)}" for i, (_, row) in enumerate(controls.iterrows())]
        fig.text(0.835, 0.5, "Controls\n\n" + "\n".join(lines), va="center", ha="left",
                 fontsize=9, color=INK, linespacing=1.6)
        fig.tight_layout(rect=(0, 0, 0.82, 1))
    else:
        fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--target-dir", type=Path, required=True,
                        help="Target directory, e.g. <TASK_ROOT>/<PROTEIN> (contains fine_screening/)")
    parser.add_argument("--target-name", required=True, help="Target name used in the plot filename")
    parser.add_argument("--scores", type=Path, default=None,
                        help="Scores CSV (default: <target-dir>/fine_screening/compiled_scores.csv)")
    parser.add_argument("--controls", type=Path, default=None,
                        help="Controls CSV with SMILES + the same metric columns (+ optional name)")
    parser.add_argument("--output-dir", type=Path, default=None,
                        help="Output directory (default: <target-dir>/fine_screening)")
    args = parser.parse_args()

    target_dir = args.target_dir.expanduser()
    scores_path = args.scores or target_dir / "fine_screening" / "compiled_scores.csv"
    out_dir = args.output_dir or target_dir / "fine_screening"
    if not scores_path.is_file() or scores_path.stat().st_size == 0:
        sys.exit(f"ERROR: scores file not found or empty: {scores_path}")
    out_dir.mkdir(parents=True, exist_ok=True)

    screen = load_metrics(scores_path)
    metrics = [m for m in METRICS if m[0] in screen.columns]
    for m, _, _ in METRICS:
        if m not in screen.columns:
            print(f"WARNING: {m} not in {scores_path.name} (method not run?); skipping")
    if not metrics:
        sys.exit("ERROR: none of the summary metrics are present")

    out_csv = out_dir / "metric_summary.csv"
    add_percentiles(screen, screen, metrics).to_csv(out_csv, index=False)
    print(f"{'metric':<24}{'n':>7}{'missing':>9}{'median':>9}{'p10':>9}{'p90':>9}{'min':>9}{'max':>9}")
    for m, _, _ in metrics:
        s = screen[m]
        print(f"{m:<24}{s.notna().sum():>7}{s.isna().sum():>9}{s.median():>9.2f}"
              f"{s.quantile(0.1):>9.2f}{s.quantile(0.9):>9.2f}{s.min():>9.2f}{s.max():>9.2f}")
    print(f"Saved: {out_csv}")

    controls = None
    if args.controls:
        controls = load_metrics(args.controls.expanduser())
        for m, _, _ in metrics:
            if m not in controls.columns:
                print(f"WARNING: controls file has no {m}; leaving it blank for controls")
                controls[m] = np.nan
        controls = add_percentiles(controls, screen, metrics)
        ctrl_csv = out_dir / "metric_summary_controls.csv"
        controls.to_csv(ctrl_csv, index=False)
        print("\ncontrols: value (% of screen beaten)")
        for i, (_, row) in enumerate(controls.iterrows()):
            cells = [f"{m}={row[m]:.2f} ({row[m + '_pct']:.0f}%)" if pd.notna(row[m]) else f"{m}=NA"
                     for m, _, _ in metrics]
            print(f"  {i + 1} {control_label(row, i)}: " + ", ".join(cells))
        print(f"Saved: {ctrl_csv}")

    out_png = out_dir / f"{args.target_name}_metric_kde.png"
    if plot(screen, controls, metrics, out_png):
        print(f"Saved: {out_png}")


if __name__ == "__main__":
    main()
