#!/usr/bin/env python3
"""Bootstrap per-protein selected.csv from a seed compounds CSV.

This script creates a canonical selected.csv with a single required column:
`SMILES`. It can also emit a snapshot copy (selected_seed.csv) used by the
AF3-first pipeline for traceability.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def detect_smiles_column(df: pd.DataFrame, requested: str) -> str:
    """Return a valid smiles column name from the dataframe."""
    if requested in df.columns:
        return requested

    lowered = {col.lower(): col for col in df.columns}
    for candidate in ("smiles", "ligand_description"):
        if candidate in lowered:
            return lowered[candidate]

    raise ValueError(
        "No SMILES-like column found. Checked requested column "
        f"'{requested}' and fallbacks: smiles, ligand_description."
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Create selected.csv from a seed compounds CSV.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--seed-csv", required=True, help="Input seed compounds CSV")
    parser.add_argument(
        "--output-selected",
        required=True,
        help="Output selected.csv path (used by downstream workflows)",
    )
    parser.add_argument(
        "--output-seed-copy",
        default="",
        help="Optional path to write selected_seed.csv snapshot",
    )
    parser.add_argument(
        "--smiles-col",
        default="SMILES",
        help="Preferred SMILES column name in seed CSV",
    )
    args = parser.parse_args()

    seed_path = Path(args.seed_csv)
    selected_path = Path(args.output_selected)
    selected_path.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(seed_path)
    smiles_col = detect_smiles_column(df, args.smiles_col)

    # Keep deterministic order from the original file, remove blanks/duplicates.
    out = pd.DataFrame({"SMILES": df[smiles_col]})
    out["SMILES"] = out["SMILES"].astype("string").str.strip()
    out = out[out["SMILES"].notna() & (out["SMILES"] != "")]
    out = out.drop_duplicates(subset=["SMILES"], keep="first").reset_index(drop=True)

    if out.empty:
        raise ValueError(f"No valid SMILES rows found in {seed_path}")

    out.to_csv(selected_path, index=False)
    print(f"Wrote {len(out)} rows to {selected_path}")

    if args.output_seed_copy:
        seed_copy_path = Path(args.output_seed_copy)
        seed_copy_path.parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(seed_copy_path, index=False)
        print(f"Wrote seed snapshot to {seed_copy_path}")


if __name__ == "__main__":
    main()
