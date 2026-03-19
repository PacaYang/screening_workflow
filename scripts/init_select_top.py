#!/usr/bin/env python3

###################################
# To select the top candidates from Initial screening for the fine screening
# Updated to use Polars and read multiple files per method from directories
###################################

import argparse
import polars as pl
from pathlib import Path

def load_csvs_from_dir(dir_path, id_col, score_col, rename_score_to):
    """Load all CSVs from a directory and combine them."""
    dir_path = Path(dir_path)
    if not dir_path.exists():
        raise ValueError(f"Directory does not exist: {dir_path}")

    csv_files = [f for f in dir_path.glob("prediction_*.csv") if not f.name.endswith("failed.csv")]
    if not csv_files:
        raise ValueError(f"No CSV files found in {dir_path}")
    
    dfs = []
    for csv_file in csv_files:
        df = pl.read_csv(csv_file)
        if id_col not in df.columns:
            raise ValueError(f"{csv_file} missing id_col '{id_col}'")
        if score_col not in df.columns:
            raise ValueError(f"{csv_file} missing score_col '{score_col}'")
        if df.is_empty():
            continue
        # Cast score column to Float64 (empty CSVs read as Utf8)
        df = df.with_columns(pl.col(score_col).cast(pl.Float64))
        dfs.append(df.select([id_col, score_col]))
    
    # Concatenate all dataframes and remove duplicates
    if not dfs:
        raise ValueError(f"No non-empty CSV files found in {dir_path}")
    combined = pl.concat(dfs)
    combined = combined.rename({score_col: rename_score_to})
    combined = combined.unique(subset=[id_col])
    print(combined)

    return combined

def even_split_select(df, left_col, right_col, n, already_selected):
    """Select n rows evenly from left_col and right_col (highest scores first),
    excluding already_selected ids. Rebalance if one side is short."""
    
    # Filter out already selected
    mask_pool = ~df["__id__"].is_in(already_selected)
    
    # Left half
    n_left = n // 2
    left_df = df.filter(mask_pool & df[left_col].is_not_null()).sort(left_col, descending=True)
    chosen_left = left_df.head(n_left)["__id__"].to_list()
    
    # Right half
    mask_pool2 = mask_pool & (~df["__id__"].is_in(chosen_left))
    n_right = n - len(chosen_left)
    right_df = df.filter(mask_pool2 & df[right_col].is_not_null()).sort(right_col, descending=True)
    chosen_right = right_df.head(n_right)["__id__"].to_list()
    
    chosen = chosen_left + chosen_right
    
    # Rebalance if we still came up short (take remaining best across both)
    short = n - len(chosen)
    if short > 0:
        mask_pool3 = mask_pool & (~df["__id__"].is_in(chosen))
        both_rank = df.filter(mask_pool3).with_columns(
            pl.max_horizontal([left_col, right_col]).alias("best")
        ).sort("best", descending=True)
        chosen += both_rank.head(short)["__id__"].to_list()
    
    return chosen

def main():
    ap = argparse.ArgumentParser(
        description="Select compounds for fine screening from multiple CSVs per method.", 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    ap.add_argument("--graphdta-dir", required=True, help="Directory containing GraphDTA CSV files")
    ap.add_argument("--hmsa-dir", required=True, help="Directory containing HMSA CSV files")
    ap.add_argument("--colddta-dir", required=True, help="Directory containing ColdDTA CSV files")
    ap.add_argument("--evidti-dir", required=False, help="Directory containing EviDTI CSV files")
    ap.add_argument("--druglamp-dir", required=False, help="Directory containing DrugLAMP CSV files")
    ap.add_argument("--conplex-dir", required=False, help="Directory containing ConPLex CSV files")
    ap.add_argument("--id-col", default="SMILES", help="Identifier column")
    ap.add_argument("--graphdta-score", default="predicted_affinity", help="GraphDTA score column")
    ap.add_argument("--hmsa-score", default="label", help="HMSA score column")
    ap.add_argument("--colddta-score", default="predicted_affinity", help="ColdDTA score column")
    ap.add_argument("--evidti-score", default="evidti_score", help="EviDTI score column")
    ap.add_argument("--druglamp-score", default="druglamp_score", help="DrugLAMP score column")
    ap.add_argument("--conplex-score", default="conplex_score", help="ConPLex score column")
    ap.add_argument("--hmsa-threshold", type=float, default=0.5, help="Keep all HMSA >= this")
    ap.add_argument("--target-n", type=int, default=10000, help="Total selected size")
    ap.add_argument("--summary", required=True, help="Path to write summary.csv")
    ap.add_argument("--selected", required=True, help="Path to write selected.csv; This file will be used to map the SMILES for Boltz2 output as well.")
    args = ap.parse_args()

    # Load & standardize from directories
    g = load_csvs_from_dir(args.graphdta_dir, args.id_col, args.graphdta_score, "graphdta_score")
    h = load_csvs_from_dir(args.hmsa_dir, args.id_col, args.hmsa_score, "hmsa_score")
    c = load_csvs_from_dir(args.colddta_dir, args.id_col, args.colddta_score, "colddta_score")

    # Load new methods (if provided)
    e = None
    d = None
    cp = None

    if args.evidti_dir:
        try:
            e = load_csvs_from_dir(args.evidti_dir, args.id_col, args.evidti_score, "evidti_score")
            print(f"Loaded EviDTI data: {len(e)} compounds")
        except Exception as ex:
            print(f"Warning: Could not load EviDTI data: {ex}")

    if args.druglamp_dir:
        try:
            d = load_csvs_from_dir(args.druglamp_dir, args.id_col, args.druglamp_score, "druglamp_score")
            print(f"Loaded DrugLAMP data: {len(d)} compounds")
            if len(d) == 0:
                print("Warning: DrugLAMP returned 0 compounds, skipping")
                d = None
        except Exception as ex:
            print(f"Warning: Could not load DrugLAMP data: {ex}")

    if args.conplex_dir:
        try:
            cp = load_csvs_from_dir(args.conplex_dir, args.id_col, args.conplex_score, "conplex_score")
            print(f"Loaded ConPLex data: {len(cp)} compounds")
            if len(cp) == 0:
                print("Warning: ConPLex returned 0 compounds, skipping")
                cp = None
        except Exception as ex:
            print(f"Warning: Could not load ConPLex data: {ex}")

    # Outer merge to keep all rows we have scores for
    df = g.join(h, on=args.id_col, how="full", coalesce=True).join(c, on=args.id_col, how="full", coalesce=True)

    # Add new methods if available
    if e is not None:
        df = df.join(e, on=args.id_col, how="full", coalesce=True)
    if d is not None:
        df = df.join(d, on=args.id_col, how="full", coalesce=True)
    if cp is not None:
        df = df.join(cp, on=args.id_col, how="full", coalesce=True)

    # Normalize id column name internally
    df = df.rename({args.id_col: "__id__"})

    # 1) Seed with HMSA >= threshold
    h_keep = df.filter(df["hmsa_score"] >= args.hmsa_threshold)
    h_keep_ids = set(h_keep["__id__"].to_list())

    # 2) If we already exceed target, cap by highest HMSA
    if len(h_keep_ids) >= args.target_n:
        capped = (h_keep
                  .sort("hmsa_score", descending=True)
                  .head(args.target_n)
                  .select(["__id__", "hmsa_score"])
                  .with_columns(pl.lit("HMSA>=thr").alias("source")))
        
        # Write outputs
        out_summary = df.join(capped.select(["__id__", "source"]), on="__id__", how="left")
        out_summary = out_summary.rename({"__id__": args.id_col})
        
        Path(args.summary).parent.mkdir(parents=True, exist_ok=True)
        out_summary.write_csv(args.summary)
        capped.rename({"__id__": args.id_col}).write_csv(args.selected)
        return

    # 3) Fill remaining from all available methods (GraphDTA, ColdDTA, DrugLAMP, ConPLex)
    remaining_n = args.target_n - len(h_keep_ids)

    # Build list of available score columns
    score_cols = ["graphdta_score", "colddta_score"]
    if d is not None:
        score_cols.append("druglamp_score")
    if cp is not None:
        score_cols.append("conplex_score")

    # Select based on best score across all available methods
    mask_pool = ~df["__id__"].is_in(h_keep_ids)
    candidates = df.filter(mask_pool).with_columns(
        pl.max_horizontal(score_cols).alias("best_score")
    ).sort("best_score", descending=True)

    chosen_even = candidates.head(remaining_n)["__id__"].to_list()
    chosen_ids = list(h_keep_ids) + chosen_even

    # Build selected table with provenance
    sel = df.filter(df["__id__"].is_in(chosen_ids))

    # Build method list for source label
    method_names = ["GraphDTA", "ColdDTA"]
    if d is not None:
        method_names.append("DrugLAMP")
    if cp is not None:
        method_names.append("ConPLex")
    fill_label = f"FILL({'/'.join(method_names)})"

    sel = sel.with_columns(
        pl.when(pl.col("__id__").is_in(list(h_keep_ids)))
          .then(pl.lit("HMSA>=thr"))
          .otherwise(pl.lit(fill_label))
          .alias("source")
    )

    # Order by (HMSA first, then best of others)
    all_score_cols = score_cols + ["hmsa_score"]
    sel = sel.with_columns([
        (pl.col("source") == "HMSA>=thr").cast(pl.Int32).alias("_sort_h"),
        pl.max_horizontal(all_score_cols).alias("_sort_best")
    ]).sort(["_sort_h", "_sort_best"], descending=[True, True])

    # Outputs
    out_summary = df.join(sel.select(["__id__", "source"]), on="__id__", how="left")
    out_summary = out_summary.rename({"__id__": args.id_col})
    sel_out = sel.rename({"__id__": args.id_col})

    # Build output columns dynamically based on available methods
    output_cols = [args.id_col, "source", "hmsa_score", "graphdta_score", "colddta_score"]
    if d is not None:
        output_cols.append("druglamp_score")
    if cp is not None:
        output_cols.append("conplex_score")

    Path(args.summary).parent.mkdir(parents=True, exist_ok=True)
    out_summary.write_csv(args.summary)
    sel_out.select(output_cols).write_csv(args.selected)

if __name__ == "__main__":
    main()
