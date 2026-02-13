#!/usr/bin/env python3

###################################
# To select the top candidates from Initial screening for the fine screening
###################################

import argparse
import pandas as pd
from pathlib import Path

def load_csv(path, id_col, score_col, rename_score_to):
    df = pd.read_csv(path)
    if id_col not in df.columns:
        raise ValueError(f"{path} missing id_col '{id_col}'")
    if score_col not in df.columns:
        raise ValueError(f"{path} missing score_col '{score_col}'")
    df = df[[id_col, score_col]].rename(columns={score_col: rename_score_to})
    return df.drop_duplicates(subset=[id_col])

def even_split_select(df, left_col, right_col, n, already_selected):
    """Select n rows evenly from left_col and right_col (highest scores first),
    excluding already_selected ids. Rebalance if one side is short."""
    ids = df["__id__"]
    mask_pool = ~ids.isin(already_selected)

    # Left half
    n_left = n // 2
    left_df = df[mask_pool & df[left_col].notna()].sort_values(left_col, ascending=False)
    chosen_left = left_df.head(n_left)["__id__"].tolist()

    # Right half
    mask_pool2 = mask_pool & (~ids.isin(chosen_left))
    n_right = n - len(chosen_left)
    right_df = df[mask_pool2 & df[right_col].notna()].sort_values(right_col, ascending=False)
    chosen_right = right_df.head(n_right)["__id__"].tolist()

    chosen = chosen_left + chosen_right
    # Rebalance if we still came up short (take remaining best across both)
    short = n - len(chosen)
    if short > 0:
        mask_pool3 = mask_pool & (~ids.isin(chosen))
        both_rank = df[mask_pool3].assign(
            best=df[[left_col, right_col]].max(axis=1, skipna=True)
        ).sort_values("best", ascending=False)
        chosen += both_rank.head(short)["__id__"].tolist()
    return chosen

def main():
    ap = argparse.ArgumentParser(description="Select compounds for fine screening.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--graphdta", required=True, help="GraphDTA CSV")
    ap.add_argument("--hmsa", required=True, help="HMSA CSV")
    ap.add_argument("--colddta", required=True, help="ColdDTA CSV")
    ap.add_argument("--id-col", default="SMILES", help="Identifier column")
    ap.add_argument("--graphdta-score", default="predicted_affinity", help="GraphDTA score column")
    ap.add_argument("--hmsa-score", default="label", help="HMSA score column")
    ap.add_argument("--colddta-score", default="predicted_affinity", help="ColdDTA score column")
    ap.add_argument("--hmsa-threshold", type=float, default=0.5, help="Keep all HMSA >= this")
    ap.add_argument("--target-n", type=int, default=10000, help="Total selected size")
    ap.add_argument("--summary", required=True, help="Path to write summary.csv")
    ap.add_argument("--selected", required=True, help="Path to write selected.csv; This filed will be used to map the SMILES for Boltz2 output as well.")
    args = ap.parse_args()

    # Load & standardize
    g = load_csv(args.graphdta, args.id_col, args.graphdta_score, "graphdta_score")
    h = load_csv(args.hmsa,     args.id_col, args.hmsa_score,     "hmsa_score")
    c = load_csv(args.colddta,  args.id_col, args.colddta_score,  "colddta_score")

    # Outer merge to keep all rows we have scores for
    df = g.merge(h, on=args.id_col, how="outer").merge(c, on=args.id_col, how="outer")
    # Normalize id column name internally
    df = df.rename(columns={args.id_col: "__id__"})

    # 1) Seed with HMSA >= threshold
    h_keep = df[df["hmsa_score"].ge(args.hmsa_threshold, fill_value=False)].copy()
    h_keep_ids = set(h_keep["__id__"].tolist())

    # 2) If we already exceed target, cap by highest HMSA
    if len(h_keep_ids) >= args.target_n:
        capped = (h_keep
                  .sort_values("hmsa_score", ascending=False)
                  .head(args.target_n)[["__id__", "hmsa_score"]]
                  .assign(source="HMSA>=thr"))
        # Write outputs
        out_summary = df.merge(capped[["__id__", "source"]], on="__id__", how="left")
        out_summary = out_summary.rename(columns={"__id__": args.id_col})
        Path(args.summary).parent.mkdir(parents=True, exist_ok=True)
        out_summary.to_csv(args.summary, index=False)
        capped.rename(columns={"__id__": args.id_col}).to_csv(args.selected, index=False)
        return

    # 3) Fill remaining evenly from GraphDTA and ColdDTA
    remaining_n = args.target_n - len(h_keep_ids)
    chosen_even = even_split_select(
        df=df,
        left_col="graphdta_score",
        right_col="colddta_score",
        n=remaining_n,
        already_selected=h_keep_ids
    )
    chosen_ids = list(h_keep_ids) + chosen_even

    # Build selected table with provenance
    sel = df[df["__id__"].isin(chosen_ids)].copy()
    sel["source"] = "FILL(GraphDTA/ColdDTA)"
    sel.loc[sel["__id__"].isin(h_keep_ids), "source"] = "HMSA>=thr"

    # Order by (HMSA first, then best of others)
    sel = sel.assign(
        _sort_h=sel["source"].eq("HMSA>=thr").astype(int),
        _sort_best=sel[["graphdta_score", "colddta_score", "hmsa_score"]].max(axis=1, skipna=True)
    ).sort_values(["_sort_h", "_sort_best"], ascending=[False, False])

    # Outputs
    out_summary = df.merge(sel[["__id__", "source"]], on="__id__", how="left")
    out_summary = out_summary.rename(columns={"__id__": args.id_col})
    sel_out = sel.rename(columns={"__id__": args.id_col})

    Path(args.summary).parent.mkdir(parents=True, exist_ok=True)
    out_summary.to_csv(args.summary, index=False)
    sel_out[[args.id_col, "source", "hmsa_score", "graphdta_score", "colddta_score"]].to_csv(
        args.selected, index=False
    )

if __name__ == "__main__":
    main()
