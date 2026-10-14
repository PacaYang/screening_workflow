#!/usr/bin/env python3
"""Merge available fine-screening method summaries into one CSV."""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd


METHOD_FILES = {
    "AF3": ("AF3/summary.csv", "af3"),
    "Boltz2": ("Boltz2/summary.csv", "boltz2"),
    "Vina": ("Vina/results.csv", "vina"),
    "RoseTTAFold": ("RoseTTAFold/summary.csv", "rfaa"),
    "Nesso": ("Nesso/summary.csv", "nesso"),
    "Protenix": ("Protenix/summary.csv", "protenix"),
    "OpenFold3": ("OpenFold3/summary.csv", "openfold3"),
}
FOLDER_RE = re.compile(r"^(?:job|compound)_(\d+)$")


def read_scores_csv(path):
    frame = pd.read_csv(path)
    return frame.loc[:, ~frame.columns.astype(str).str.startswith("Unnamed:")]


def normalize_smiles(frame, path):
    if "SMILES" not in frame:
        raise ValueError(f"summary has no SMILES column: {path}")
    frame = frame.copy()
    frame["SMILES"] = frame["SMILES"].fillna("").astype(str).str.strip()
    if (frame["SMILES"] == "").any():
        raise ValueError(f"summary has blank SMILES: {path}")
    return frame


def selected_frame(path):
    selected = normalize_smiles(read_scores_csv(path), path)
    selected.insert(0, "compound_index", range(len(selected)))
    return selected


def map_folder_indices(frame, path):
    if "compound_index" in frame:
        indices = pd.to_numeric(frame["compound_index"], errors="raise")
        if indices.isna().any() or (indices % 1 != 0).any():
            raise ValueError(f"invalid compound_index in {path}")
        return indices.astype(int)
    if "folder" not in frame:
        raise ValueError(f"summary lacks compound_index/folder: {path}")
    values = frame["folder"].astype(str).str.extract(FOLDER_RE, expand=False)
    if values.isna().any():
        raise ValueError(f"cannot map folder to compound index in {path}")
    return values.astype(int)


def map_vina_indices(frame, target_dir, path):
    input_dir = target_dir / "fine_screening" / "Vina" / "input"
    files = sorted(input_dir.glob("input_part_*.csv"),
                   key=lambda p: int(p.stem.rsplit("_", 1)[1]))
    partitioned = bool(files)
    if not files:
        files = sorted(input_dir.glob("*.csv"))
    if frame.empty:
        return pd.Series(dtype=int, index=frame.index)
    if not files:
        raise ValueError(f"Vina input CSV files not found under {input_dir}")
    lookup = {}
    offset = 0
    for fallback_part, file in enumerate(files):
        part = int(file.stem.rsplit("_", 1)[1]) if partitioned else fallback_part
        table = pd.read_csv(file)
        for local_index in range(len(table)):
            lookup[(part, local_index)] = offset + local_index
        offset += len(table)
    try:
        indices = [lookup[(int(folder), int(idx))]
                   for folder, idx in zip(frame["folder"], frame["idx"])]
    except (KeyError, ValueError, TypeError) as error:
        raise ValueError(f"cannot map Vina folder/idx in {path}: {error}") from error
    return pd.Series(indices, index=frame.index, dtype=int)


def prepare_method(frame, method, prefix, target_dir, path, selected=None):
    frame = normalize_smiles(frame, path)
    if method == "Vina":
        frame["affinity"] = pd.to_numeric(frame["affinity"], errors="coerce")
        frame = frame.dropna(subset=["affinity"]).copy()
    if selected is not None:
        if method == "Vina":
            indices = map_vina_indices(frame, target_dir, path)
        else:
            indices = map_folder_indices(frame, path)
        frame = frame.copy()
        frame["compound_index"] = indices
        if frame["compound_index"].duplicated().any():
            if method != "Vina":
                raise ValueError(f"duplicate compound_index in {path}")
            frame = frame.sort_values("affinity", ascending=True).drop_duplicates("compound_index")
        if not frame["compound_index"].isin(selected["compound_index"]).all():
            raise ValueError(f"compound index outside selected.csv in {path}")
        checked = frame.merge(selected[["compound_index", "SMILES"]],
                              on="compound_index", suffixes=("", "_selected"), validate="one_to_one")
        if (checked["SMILES"] != checked["SMILES_selected"]).any():
            raise ValueError(f"SMILES mismatch between {path} and selected.csv")
        frame = frame.drop(columns=["SMILES"])
        key = "compound_index"
    else:
        if method == "Vina":
            frame = frame.sort_values("affinity", ascending=True).drop_duplicates("SMILES")
        if frame["SMILES"].duplicated().any():
            raise ValueError(f"duplicate SMILES in {path}; selected.csv is required for index joins")
        key = "SMILES"
    return frame.rename(columns={col: f"{prefix}_{col}" for col in frame if col != key})


def compile_scores(target_dir, rfaa_path=None, skip_methods=()):
    fine_dir = target_dir / "fine_screening"
    selected_path = target_dir / "initial_screening" / "selected.csv"
    selected = selected_frame(selected_path) if selected_path.is_file() else None
    merged = selected.copy() if selected is not None else None
    loaded = []
    for method, (relative, prefix) in METHOD_FILES.items():
        if method in skip_methods:
            continue
        path = Path(rfaa_path) if method == "RoseTTAFold" and rfaa_path else fine_dir / relative
        if not path.is_file() or path.stat().st_size == 0:
            continue
        frame = read_scores_csv(path)
        prepared = prepare_method(frame, method, prefix, target_dir, path, selected)
        if merged is None:
            merged = prepared
        else:
            key = "compound_index" if selected is not None else "SMILES"
            merged = merged.merge(prepared, on=key, how="left" if selected is not None else "outer",
                                  validate="one_to_one")
        loaded.append(method)
    if not loaded:
        raise ValueError(f"no method summaries found under {fine_dir}")
    if selected is None:
        merged = merged[["SMILES"] + [col for col in merged if col != "SMILES"]]
    return merged, loaded


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-dir", type=Path, required=True)
    parser.add_argument("--target-name", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--rfaa-path", type=Path, default=None)
    parser.add_argument("--skip-method", action="append", choices=METHOD_FILES, default=[])
    args = parser.parse_args()
    try:
        frame, methods = compile_scores(args.target_dir.expanduser(), args.rfaa_path, args.skip_method)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(args.output, index=False)
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    print(f"Wrote {len(frame)} rows for {args.target_name} to {args.output}; methods: {', '.join(methods)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
