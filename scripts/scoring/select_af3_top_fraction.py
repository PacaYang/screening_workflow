#!/usr/bin/env python3
"""Select top AF3 compounds using pLDDT and chain_pair_pae_min criteria.

Selection policy:
1) Keep rows with plddt >= threshold.
2) Ignore rows where chain_pair_pae_min == -1.
3) Rank by chain_pair_pae_min ascending (lower is better),
   tie-break by plddt descending.
4) Select max(1, ceil(top_fraction * eligible_count)).
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Select top fraction of AF3 results for downstream screening.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--af3-summary", required=True, help="AF3 summary.csv path")
    parser.add_argument(
        "--selected-output",
        required=True,
        help="Output selected.csv path (used by downstream workflows)",
    )
    parser.add_argument(
        "--top-output",
        default="",
        help="Optional detailed output path for selected AF3 rows",
    )
    parser.add_argument(
        "--plddt-threshold",
        type=float,
        default=70.0,
        help="Minimum pLDDT filter",
    )
    parser.add_argument(
        "--top-fraction",
        type=float,
        default=0.30,
        help="Top fraction to keep from eligible compounds",
    )
    args = parser.parse_args()

    if args.top_fraction <= 0 or args.top_fraction > 1:
        raise ValueError("--top-fraction must be in (0, 1].")

    summary_path = Path(args.af3_summary)
    selected_output = Path(args.selected_output)
    selected_output.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(summary_path)

    required_cols = {"SMILES", "plddt", "chain_pair_pae_min"}
    missing = required_cols.difference(df.columns)
    if missing:
        raise ValueError(
            f"Missing required AF3 columns in {summary_path}: {sorted(missing)}"
        )

    work = df.copy()
    work["SMILES"] = work["SMILES"].astype("string").str.strip()
    work["plddt"] = pd.to_numeric(work["plddt"], errors="coerce")
    work["chain_pair_pae_min"] = pd.to_numeric(
        work["chain_pair_pae_min"], errors="coerce"
    )

    eligible = work[
        (work["SMILES"].notna())
        & (work["SMILES"] != "")
        & (work["plddt"] >= args.plddt_threshold)
        & (work["chain_pair_pae_min"].notna())
        & (work["chain_pair_pae_min"] != -1)
    ].copy()

    if eligible.empty:
        raise ValueError(
            "No eligible AF3 rows after filtering. "
            f"Threshold={args.plddt_threshold}, invalid chain_pair_pae_min removed."
        )

    eligible = eligible.sort_values(
        by=["chain_pair_pae_min", "plddt"], ascending=[True, False]
    )

    # One row per SMILES for downstream selected.csv.
    ranked = eligible.drop_duplicates(subset=["SMILES"], keep="first").reset_index(
        drop=True
    )

    n_eligible = len(ranked)
    n_select = max(1, int(math.ceil(args.top_fraction * n_eligible)))
    selected = ranked.head(n_select).copy()
    selected["selection_rank"] = range(1, len(selected) + 1)
    selected["source"] = "AF3_top_fraction"

    # Compatibility file consumed by downstream scripts.
    selected[["SMILES"]].to_csv(selected_output, index=False)
    print(
        f"Selected {len(selected)}/{n_eligible} compounds "
        f"({args.top_fraction:.2%}) -> {selected_output}"
    )

    if args.top_output:
        top_output = Path(args.top_output)
        top_output.parent.mkdir(parents=True, exist_ok=True)
        cols = ["SMILES", "selection_rank", "chain_pair_pae_min", "plddt", "source"]
        selected[cols].rename(
            columns={
                "chain_pair_pae_min": "af3_chain_pair_pae_min",
                "plddt": "af3_plddt",
            }
        ).to_csv(top_output, index=False)
        print(f"Wrote detailed selection to {top_output}")


if __name__ == "__main__":
    main()
