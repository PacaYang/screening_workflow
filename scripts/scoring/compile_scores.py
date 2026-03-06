#!/usr/bin/env python3
"""
Compile screening scores from AF3, Boltz2, and Vina into a single CSV file.

This script merges scores from multiple structure prediction and docking methods,
using SMILES as the primary key. For Vina, it selects the best (most negative)
affinity across all docking boxes.
"""

import argparse
import pandas as pd
from pathlib import Path
import sys


def load_af3_scores(af3_path):
    """Load AF3 scores and add method prefix to columns."""
    print(f"Loading AF3 scores from {af3_path}")
    df = pd.read_csv(af3_path, index_col=0)

    # Rename columns with af3_ prefix (except SMILES)
    columns_to_rename = {col: f'af3_{col}' for col in df.columns if col != 'SMILES'}
    df = df.rename(columns=columns_to_rename)

    print(f"  Loaded {len(df)} AF3 entries")
    return df


def load_boltz2_scores(boltz2_path):
    """Load Boltz2 scores and add method prefix to columns."""
    print(f"Loading Boltz2 scores from {boltz2_path}")
    df = pd.read_csv(boltz2_path)

    # Rename columns with boltz2_ prefix (except SMILES)
    columns_to_rename = {col: f'boltz2_{col}' for col in df.columns if col != 'SMILES'}
    df = df.rename(columns=columns_to_rename)

    print(f"  Loaded {len(df)} Boltz2 entries")
    return df


def load_vina_scores(vina_path):
    """Load Vina scores and select best affinity per SMILES."""
    print(f"Loading Vina scores from {vina_path}")
    df = pd.read_csv(vina_path, index_col=0)

    # Group by SMILES and take the row with minimum (most negative) affinity
    print(f"  Loaded {len(df)} Vina docking results")
    best_df = df.loc[df.groupby('SMILES')['affinity'].idxmin()]

    # Rename columns with vina_ prefix (except SMILES)
    columns_to_rename = {col: f'vina_{col}' for col in best_df.columns if col != 'SMILES'}
    best_df = best_df.rename(columns=columns_to_rename)

    print(f"  Selected best affinity for {len(best_df)} unique compounds")
    return best_df


def merge_scores(af3_df, boltz2_df, vina_df):
    """Merge all score dataframes on SMILES using outer join."""
    print("\nMerging scores...")

    # Start with AF3
    merged = af3_df.copy()
    print(f"  Starting with {len(merged)} AF3 entries")

    # Merge with Boltz2
    merged = merged.merge(boltz2_df, on='SMILES', how='outer')
    print(f"  After Boltz2 merge: {len(merged)} entries")

    # Merge with Vina
    merged = merged.merge(vina_df, on='SMILES', how='outer')
    print(f"  After Vina merge: {len(merged)} entries")

    # Move SMILES to first column
    cols = ['SMILES'] + [col for col in merged.columns if col != 'SMILES']
    merged = merged[cols]

    return merged


def print_summary(df, target_name):
    """Print summary statistics for the merged dataframe."""
    print(f"\n{'='*60}")
    print(f"Summary for {target_name}")
    print(f"{'='*60}")

    print(f"\nTotal unique compounds: {len(df)}")

    # Count non-null values for each method
    af3_cols = [col for col in df.columns if col.startswith('af3_')]
    boltz2_cols = [col for col in df.columns if col.startswith('boltz2_')]
    vina_cols = [col for col in df.columns if col.startswith('vina_')]

    if af3_cols:
        af3_count = df[af3_cols[0]].notna().sum()
        print(f"Compounds with AF3 scores: {af3_count}")

    if boltz2_cols:
        boltz2_count = df[boltz2_cols[0]].notna().sum()
        print(f"Compounds with Boltz2 scores: {boltz2_count}")

    if vina_cols:
        vina_count = df[vina_cols[0]].notna().sum()
        print(f"Compounds with Vina scores: {vina_count}")

    # Count compounds with all scores
    all_methods_mask = True
    if af3_cols:
        all_methods_mask &= df[af3_cols[0]].notna()
    if boltz2_cols:
        all_methods_mask &= df[boltz2_cols[0]].notna()
    if vina_cols:
        all_methods_mask &= df[vina_cols[0]].notna()

    complete_count = all_methods_mask.sum()
    print(f"Compounds with complete scores (all methods): {complete_count}")

    print(f"\nFirst 5 rows:")
    print(df.head().to_string())
    print(f"\n{'='*60}\n")


def main():
    parser = argparse.ArgumentParser(
        description='Compile screening scores from AF3, Boltz2, and Vina'
    )
    parser.add_argument(
        '--target-dir',
        required=True,
        help='Target directory containing fine_screening subdirectory'
    )
    parser.add_argument(
        '--target-name',
        required=True,
        help='Target name (e.g., IL6RA, IL6RB)'
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Output CSV file path'
    )

    args = parser.parse_args()

    # Expand paths
    target_dir = Path(args.target_dir).expanduser()
    output_path = Path(args.output).expanduser()

    # Define input paths
    af3_path = target_dir / 'fine_screening' / 'AF3' / 'summary.csv'
    boltz2_path = target_dir / 'fine_screening' / 'Boltz2' / 'summary.csv'
    vina_path = target_dir / 'fine_screening' / 'Vina' / 'results.csv'

    # Check that all input files exist
    missing_files = []
    for path in [af3_path, boltz2_path, vina_path]:
        if not path.exists():
            missing_files.append(str(path))

    if missing_files:
        print("ERROR: Missing input files:", file=sys.stderr)
        for f in missing_files:
            print(f"  {f}", file=sys.stderr)
        sys.exit(1)

    # Load scores
    try:
        af3_df = load_af3_scores(af3_path)
        boltz2_df = load_boltz2_scores(boltz2_path)
        vina_df = load_vina_scores(vina_path)
    except Exception as e:
        print(f"ERROR loading scores: {e}", file=sys.stderr)
        sys.exit(1)

    # Merge scores
    try:
        merged_df = merge_scores(af3_df, boltz2_df, vina_df)
    except Exception as e:
        print(f"ERROR merging scores: {e}", file=sys.stderr)
        sys.exit(1)

    # Save output
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        merged_df.to_csv(output_path, index=False)
        print(f"\nSaved compiled scores to: {output_path}")
    except Exception as e:
        print(f"ERROR saving output: {e}", file=sys.stderr)
        sys.exit(1)

    # Print summary
    print_summary(merged_df, args.target_name)


if __name__ == '__main__':
    main()
