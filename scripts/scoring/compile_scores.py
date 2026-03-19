#!/usr/bin/env python3
"""
Compile screening scores from AF3, Boltz2, Vina, and RoseTTAFold into a single CSV file.

This script merges scores from multiple structure prediction and docking methods,
using SMILES as the primary key. For Vina, it selects the best (most negative)
affinity across all docking boxes. RoseTTAFold is optional.
"""

import argparse
import pandas as pd
from pathlib import Path
import sys

def read_scores_csv(path):
    """Read a score CSV while normalizing legacy pandas index columns."""
    df = pd.read_csv(path)
    unnamed_cols = [c for c in df.columns if str(c).startswith('Unnamed:')]
    if unnamed_cols:
        df = df.drop(columns=unnamed_cols)
    return df


def load_af3_scores(af3_path):
    """Load AF3 scores and add method prefix to columns."""
    print(f"Loading AF3 scores from {af3_path}")
    df = read_scores_csv(af3_path)
    if 'SMILES' not in df.columns:
        raise ValueError(f"AF3 summary missing SMILES column: {af3_path}")

    # Rename columns with af3_ prefix (except SMILES)
    columns_to_rename = {col: f'af3_{col}' for col in df.columns if col != 'SMILES'}
    df = df.rename(columns=columns_to_rename)

    print(f"  Loaded {len(df)} AF3 entries")
    return df


def load_boltz2_scores(boltz2_path):
    """Load Boltz2 scores and add method prefix to columns."""
    print(f"Loading Boltz2 scores from {boltz2_path}")
    df = read_scores_csv(boltz2_path)
    if 'SMILES' not in df.columns:
        raise ValueError(f"Boltz2 summary missing SMILES column: {boltz2_path}")

    # Rename columns with boltz2_ prefix (except SMILES)
    columns_to_rename = {col: f'boltz2_{col}' for col in df.columns if col != 'SMILES'}
    df = df.rename(columns=columns_to_rename)

    print(f"  Loaded {len(df)} Boltz2 entries")
    return df


def load_rosettafold_scores(rfaa_path):
    """Load RoseTTAFold scores and add method prefix to columns."""
    print(f"Loading RoseTTAFold scores from {rfaa_path}")
    df = read_scores_csv(rfaa_path)
    if 'SMILES' not in df.columns:
        raise ValueError(f"RoseTTAFold summary missing SMILES column: {rfaa_path}")

    columns_to_rename = {col: f'rfaa_{col}' for col in df.columns if col != 'SMILES'}
    df = df.rename(columns=columns_to_rename)

    print(f"  Loaded {len(df)} RoseTTAFold entries")
    return df


def load_vina_scores(vina_path):
    """Load Vina scores and select best affinity per SMILES."""
    print(f"Loading Vina scores from {vina_path}")
    df = read_scores_csv(vina_path)
    expected_cols = ['folder', 'idx', 'box', 'affinity', 'SMILES']
    for col in expected_cols:
        if col not in df.columns:
            df[col] = pd.Series(dtype='object')

    print(f"  Loaded {len(df)} Vina docking results")
    if df.empty:
        best_df = df[expected_cols].copy()
    else:
        df['affinity'] = pd.to_numeric(df['affinity'], errors='coerce')
        # Group by SMILES and take the row with minimum (most negative) affinity
        best_df = df.loc[df.groupby('SMILES')['affinity'].idxmin()]

    # Rename columns with vina_ prefix (except SMILES)
    columns_to_rename = {col: f'vina_{col}' for col in best_df.columns if col != 'SMILES'}
    best_df = best_df.rename(columns=columns_to_rename)

    print(f"  Selected best affinity for {len(best_df)} unique compounds")
    return best_df


def merge_scores(af3_df, boltz2_df, vina_df, rfaa_df=None):
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

    # Merge with RoseTTAFold (optional)
    if rfaa_df is not None:
        merged = merged.merge(rfaa_df, on='SMILES', how='outer')
        print(f"  After RoseTTAFold merge: {len(merged)} entries")

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

    rfaa_cols = [col for col in df.columns if col.startswith('rfaa_')]
    if rfaa_cols:
        rfaa_count = df[rfaa_cols[0]].notna().sum()
        print(f"Compounds with RoseTTAFold scores: {rfaa_count}")

    # Count compounds with all scores
    all_methods_mask = pd.Series(True, index=df.index)
    if af3_cols:
        all_methods_mask &= df[af3_cols[0]].notna()
    if boltz2_cols:
        all_methods_mask &= df[boltz2_cols[0]].notna()
    if vina_cols:
        all_methods_mask &= df[vina_cols[0]].notna()
    if rfaa_cols:
        all_methods_mask &= df[rfaa_cols[0]].notna()

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
    parser.add_argument(
        '--rfaa-path',
        default=None,
        help='Path to RoseTTAFold summary.csv (optional)'
    )

    args = parser.parse_args()

    # Expand paths
    target_dir = Path(args.target_dir).expanduser()
    output_path = Path(args.output).expanduser()

    # Define input paths
    af3_path = target_dir / 'fine_screening' / 'AF3' / 'summary.csv'
    boltz2_path = target_dir / 'fine_screening' / 'Boltz2' / 'summary.csv'
    vina_path = target_dir / 'fine_screening' / 'Vina' / 'results.csv'
    rfaa_path = Path(args.rfaa_path).expanduser() if args.rfaa_path else \
                target_dir / 'fine_screening' / 'RoseTTAFold' / 'summary.csv'

    # Check that required input files exist
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
        if rfaa_path.exists() and rfaa_path.stat().st_size > 1:
            rfaa_df = load_rosettafold_scores(rfaa_path)
        else:
            rfaa_df = None
            print(f"RoseTTAFold summary not found or empty at {rfaa_path}, skipping")
    except Exception as e:
        print(f"ERROR loading scores: {e}", file=sys.stderr)
        sys.exit(1)

    # Merge scores
    try:
        merged_df = merge_scores(af3_df, boltz2_df, vina_df, rfaa_df)
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
