#!/usr/bin/env python3
"""Select AF3-prioritized compounds for downstream fine-screening.

Rules:
1) Keep rows with plddt > threshold
2) Ignore rows with chain_pair_pae_min == -1
3) Rank by chain_pair_pae_min ascending (smaller is better)
4) Keep top fraction (default 50%)

Fallback policy:
- If strict selection is empty, keep all rows with valid chain_pair_pae_min != -1.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Select AF3 top compounds",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--af3-summary", required=True, help="AF3 summary CSV path")
    parser.add_argument("--selected-all", required=True, help="Input compounds CSV used for AF3 all-compounds run")
    parser.add_argument("--output", required=True, help="Output selected compounds CSV path")
    parser.add_argument("--metrics-output", default=None, help="Optional detailed metrics CSV path")
    parser.add_argument("--smiles-col", default="SMILES", help="SMILES column name in --selected-all")
    parser.add_argument("--plddt-threshold", type=float, default=70.0, help="Strict filter threshold for pLDDT")
    parser.add_argument("--top-fraction", type=float, default=0.5, help="Fraction to keep after strict filtering")
    return parser.parse_args()


def normalize_summary(summary_path: Path) -> pd.DataFrame:
    df = pd.read_csv(summary_path)
    required = ["SMILES", "chain_pair_pae_min", "plddt"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"AF3 summary missing required columns: {missing}")

    out = df[["SMILES", "chain_pair_pae_min", "plddt"]].copy()
    out = out.dropna(subset=["SMILES"])
    out["SMILES"] = out["SMILES"].astype(str)
    out["chain_pair_pae_min"] = pd.to_numeric(out["chain_pair_pae_min"], errors="coerce")
    out["plddt"] = pd.to_numeric(out["plddt"], errors="coerce")

    # Consolidate by SMILES in case of duplicates, preferring best PAE and best pLDDT
    out = (
        out.groupby("SMILES", as_index=False)
        .agg(
            chain_pair_pae_min=("chain_pair_pae_min", "min"),
            plddt=("plddt", "max"),
        )
    )
    return out


def select_compounds(summary_df: pd.DataFrame, plddt_threshold: float, top_fraction: float) -> tuple[pd.DataFrame, str]:
    if top_fraction <= 0 or top_fraction > 1:
        raise ValueError("--top-fraction must be in (0, 1]")

    valid_pae = summary_df[
        summary_df["chain_pair_pae_min"].notna() & (summary_df["chain_pair_pae_min"] != -1)
    ].copy()

    if valid_pae.empty:
        raise ValueError("No AF3 rows have valid chain_pair_pae_min (!= -1).")

    strict = valid_pae[valid_pae["plddt"] > plddt_threshold].copy()
    strict = strict.sort_values("chain_pair_pae_min", ascending=True)

    if not strict.empty:
        keep_n = max(1, math.ceil(len(strict) * top_fraction))
        selected = strict.head(keep_n).copy()
        policy = f"strict(plddt>{plddt_threshold})_top_{top_fraction:.2f}"
        return selected, policy

    # User-selected fallback: keep all with valid PAE
    selected = valid_pae.sort_values("chain_pair_pae_min", ascending=True).copy()
    policy = "fallback_all_valid_pae"
    return selected, policy


def main() -> None:
    args = parse_args()

    summary_path = Path(args.af3_summary)
    selected_all_path = Path(args.selected_all)
    output_path = Path(args.output)
    metrics_path = Path(args.metrics_output) if args.metrics_output else output_path.with_name("af3_selection_metrics.csv")

    summary_df = normalize_summary(summary_path)
    selected_df, policy = select_compounds(summary_df, args.plddt_threshold, args.top_fraction)

    all_df = pd.read_csv(selected_all_path)
    if args.smiles_col not in all_df.columns:
        raise ValueError(f"Input compounds file missing SMILES column: {args.smiles_col}")

    if args.smiles_col != "SMILES":
        all_df = all_df.rename(columns={args.smiles_col: "SMILES"})

    all_df["SMILES"] = all_df["SMILES"].astype(str)

    # Preserve AF3 ranking order by mapping SMILES to rank
    selected_df = selected_df.reset_index(drop=True)
    rank_map = {sm: i for i, sm in enumerate(selected_df["SMILES"].tolist())}

    out_df = all_df[all_df["SMILES"].isin(rank_map.keys())].copy()
    out_df["af3_rank"] = out_df["SMILES"].map(rank_map)
    out_df = out_df.sort_values("af3_rank", ascending=True).drop(columns=["af3_rank"])

    if out_df.empty:
        raise ValueError("AF3 selection produced zero rows after joining to selected-all compounds.")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    metrics_path.parent.mkdir(parents=True, exist_ok=True)

    out_df.to_csv(output_path, index=False)

    metrics = summary_df.copy()
    metrics["selected"] = metrics["SMILES"].isin(set(out_df["SMILES"]))
    metrics["selection_policy"] = policy
    metrics.to_csv(metrics_path, index=False)

    print(f"Selection policy: {policy}")
    print(f"AF3 summary rows: {len(summary_df)}")
    print(f"Selected rows: {len(out_df)}")
    print(f"Wrote selected CSV: {output_path}")
    print(f"Wrote metrics CSV: {metrics_path}")


if __name__ == "__main__":
    main()
