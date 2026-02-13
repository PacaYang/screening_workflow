#!/usr/bin/env python3
"""
make_input_csv.py

Combine a single protein sequence with a list of ligands (SMILES)
to produce multiple input_i.csv for screening methods.

Usage:
    python make_input_csv.py \
        --sequences /shared/task_3/Input/sequences.csv \
        --smiles /shared/task_3/Input/compounds_smiles.csv \
        --protein protein1 \
        --outdir /shared/task_3/protein1/initial_screen/inputs/  
"""

import argparse
import pandas as pd
import sys

def main():
    parser = argparse.ArgumentParser(description="Generate per-protein input.csv")
    parser.add_argument("--sequences", required=True, help="CSV with columns: name,sequence")
    parser.add_argument("--smiles", required=True, help="CSV with compound SMILES")
    parser.add_argument("--smi-col", default="SMILES", help="col name of SMILES")
    parser.add_argument("--protein", required=True, help="Protein name to select from sequence file")
    parser.add_argument("--outdir", required=True, help="Output folder for all csv files")
    parser.add_argument("--chunk", default=100000, help="The splitted csv file size")
    args = parser.parse_args()
  
    # --- Read sequences ---
    try:
        seq_df = pd.read_csv(args.sequences)
    except Exception as e:
        sys.exit(f"Error reading sequences file {args.sequences}: {e}")

    if "name" not in seq_df.columns or "sequence" not in seq_df.columns:
        sys.exit("Sequences file must have columns: name, sequence")

    if args.protein not in seq_df["name"].values:
        sys.exit(f"Protein '{args.protein}' not found in {args.sequences}")

    seq = seq_df.loc[seq_df["name"] == args.protein, "sequence"].iloc[0]

    # --- Read SMILES ---
    try:
        smiles_df = pd.read_csv(args.smiles)
    except Exception as e:
        sys.exit(f"Error reading SMILES file {args.smiles}: {e}")

    if args.smi_col not in smiles_df.columns:
        sys.exit(f"SMILES file must have a {args.smi_col} column, or change --smi-col option")

    # --- Build combined dataframe ---
    # replicate the protein sequence for every ligand
    out_df = smiles_df[[args.smi_col]]
    out_df.insert(1, "sequence", seq)
    out_df.insert(1, "label", 0)
    
    # --- Save ---
    for i in range(0, len(out_df), args.chunk):
        chunk = out_df.iloc[i:i + args.chunk]
        chunk.to_csv(f"{args.outdir}/input_{i}.csv", index=False)
    
    print(f"[OK] Wrote {len(out_df)} rows to {args.outdir}")

if __name__ == "__main__":
    main()
