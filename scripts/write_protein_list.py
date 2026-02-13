#!/usr/bin/env python3
"""
write_protein_list.py

Discover protein names from a sequences CSV, create per-protein folders,
and write an index JSON consumed by Snakemake for dynamic expansion.

Example:
    python write_protein_list.py \
        --sequences /shared/task_3/Input/sequences.csv \
        --task-root /shared/task_3 \
        --out-subdir initial_screen \
        --emit-index /shared/task_3/.inputs_index.json
"""

import argparse
import json
import os
import sys
from collections import OrderedDict

import pandas as pd


def fail(msg: str, code: int = 1):
    print(f"[ERROR] {msg}", file=sys.stderr)
    sys.exit(code)


def warn(msg: str):
    print(f"[WARN] {msg}", file=sys.stderr)


def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)


def is_path_like_ok(segment: str) -> bool:
    # Conservative check for common filesystem-unsafe chars
    bad = set('/\0')
    return not any(c in bad for c in segment)


def main():
    ap = argparse.ArgumentParser(description="Write protein list index and scaffold folders.")
    ap.add_argument("--sequences", required=True, help="CSV with columns: name,sequence")
    ap.add_argument("--task-root", required=True, help="Root directory of the task")
    ap.add_argument("--out-subdir", required=True, help="Per-protein subdir to create (e.g., initial_screen)")
    ap.add_argument("--emit-index", required=True, help="Output JSON path, e.g., <task_root>/.inputs_index.json")
    args = ap.parse_args()

    # Validate inputs exist
    if not os.path.isfile(args.sequences):
        fail(f"Sequences file not found: {args.sequences}")


    try:
        seq_df = pd.read_csv(args.sequences)
    except Exception as e:
        fail(f"Failed to read sequences CSV: {e}")

    # Column checks
    required_cols = {"name", "sequence"}
    missing = required_cols - set(seq_df.columns)
    if missing:
        fail(f"Sequences CSV missing columns: {', '.join(sorted(missing))}")

    # Normalize and filter rows
    df = seq_df.copy()
    # Strip whitespace from names; drop completely empty names/sequences
    df["name"] = df["name"].astype(str).str.strip()
    df["sequence"] = df["sequence"].astype(str).str.strip()
    df = df[(df["name"] != "") & (df["sequence"] != "")]
    if df.empty:
        fail("No valid (name, sequence) rows found after cleaning.")

    # Preserve first occurrence order of names
    ordered_names = list(OrderedDict.fromkeys(df["name"].tolist()))

    # Warn about duplicates (same name appears multiple times)
    if df["name"].duplicated().any():
        dups = df["name"][df["name"].duplicated()].unique().tolist()
        warn(f"Duplicate protein names detected (keeping first occurrence): {dups}")

    # Sanity check path segments
    for name in ordered_names:
        if not is_path_like_ok(name):
            warn(f"Protein name may be unsafe for folder name: {name!r}. "
                 f"Consider renaming or mapping to a safe identifier.")

    # Create folder scaffolding: <task_root>/<protein>/<out-subdir>
    for name in ordered_names:
        pdir = os.path.join(args.task_root, name, args.out_subdir)
        ensure_dir(pdir)

    # Ensure index directory exists
    ensure_dir(os.path.dirname(os.path.abspath(args.emit_index)))

    # Write index JSON
    index_obj = {"proteins": ordered_names}
    try:
        with open(args.emit_index, "w", encoding="utf-8") as f:
            json.dump(index_obj, f, ensure_ascii=False)
    except Exception as e:
        fail(f"Failed to write index JSON to {args.emit_index}: {e}")

    # Friendly summary
    print(f"[OK] Wrote protein index with {len(ordered_names)} entries → {args.emit_index}")
    print(f"[OK] Ensured per-protein folders under: {args.task_root}/<protein>/{args.out_subdir}")


if __name__ == "__main__":
    main()
