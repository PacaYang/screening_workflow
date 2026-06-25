"""
Prepare receptor PDB and docking-box entries for Vina docking on the
/shared/cuteness_woundcare proteins (excluding JAK1/2/3 and TYK2).

For each non-excluded protein:
  1. Locate the first compound's AF3 _model.cif inside batch_0.tar.gz, extract
     it to a scratch dir, and write chain A to Input/protein_file/<protein>/<protein>.pdb.
  2. Read fine_screening/AF3/summary.csv binding_residues, union residue numbers
     across compounds, compute the bounding box from chain-A atoms in the new
     PDB, pad by 5 A on each side.
  3. Update Input/sequences.csv with PDB_ID (empty) and "docking box" columns,
     formatted as the nested list [[cx, cy, cz, sx, sy, sz]] expected by
     run_vina_batch.sh.
"""

import argparse
import ast
import json
import os
import re
import shutil
import sys
import tarfile
import tempfile

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scoring'))
from binding_site_utils import parse_atom_site  # noqa: E402


EXCLUDE = {"JAK1", "JAK2", "JAK3", "TYK2", "JAK_input", "Input", "logs"}
PADDING_ANG = 5.0


def discover_proteins(task_root):
    out = []
    for name in sorted(os.listdir(task_root)):
        if name in EXCLUDE:
            continue
        path = os.path.join(task_root, name)
        if not os.path.isdir(path):
            continue
        if os.path.isdir(os.path.join(path, "fine_screening", "AF3", "output")):
            out.append(name)
    return out


def find_first_cif_in_archive(tgz_path, scratch):
    """Extract the first <compound>_model.cif found in tgz_path to scratch.
    Returns the absolute path of the extracted cif, or None."""
    with tarfile.open(tgz_path, "r:gz") as tf:
        cif_member = None
        for m in tf.getmembers():
            if m.name.endswith("_model.cif"):
                cif_member = m
                break
        if cif_member is None:
            return None
        tf.extract(cif_member, path=scratch)
        return os.path.join(scratch, cif_member.name)


def cif_chain_a_to_pdb(cif_path, pdb_path):
    """Write chain A ATOM records from cif_path to pdb_path.
    Returns the list of (res_seq, x, y, z) tuples used."""
    atoms = parse_atom_site(cif_path)
    if not atoms:
        raise RuntimeError(f"No _atom_site rows parsed from {cif_path}")

    serial = 0
    coords = []
    with open(pdb_path, "w") as f:
        for a in atoms:
            if a.get("label_asym_id") != "A":
                continue
            if a.get("group_PDB", "ATOM") != "ATOM":
                continue
            try:
                x = float(a["Cartn_x"]); y = float(a["Cartn_y"]); z = float(a["Cartn_z"])
                seq = int(a["label_seq_id"])
            except (KeyError, ValueError):
                continue
            atom_name = a.get("label_atom_id", "").strip().strip('"')
            res_name = a.get("label_comp_id", "UNK")[:3]
            element = a.get("type_symbol", atom_name[:1] if atom_name else "C")[:2]
            occupancy = float(a.get("occupancy", 1.00) or 1.00)
            bfactor = float(a.get("B_iso_or_equiv", 0.00) or 0.00)
            serial += 1
            if len(atom_name) < 4:
                atom_field = f" {atom_name:<3}"
            else:
                atom_field = atom_name[:4]
            line = (
                f"ATOM  {serial:5d} {atom_field} {res_name:>3} A{seq:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{bfactor:6.2f}"
                f"          {element:>2}\n"
            )
            f.write(line)
            coords.append((seq, x, y, z))
        f.write("END\n")
    if not coords:
        raise RuntimeError(f"No chain-A ATOM records written from {cif_path}")
    return coords


def union_residue_numbers(summary_csv):
    """Return the set of residue numbers across all rows' binding_residues."""
    df = pd.read_csv(summary_csv)
    if "binding_residues" not in df.columns:
        return set()
    nums = set()
    for cell in df["binding_residues"].dropna():
        if not isinstance(cell, str) or not cell.strip():
            continue
        for token in cell.split(","):
            token = token.strip()
            if not token:
                continue
            head = token.split(":", 1)[0]
            try:
                nums.add(int(head))
            except ValueError:
                continue
    return nums


def compute_box(coords, residue_numbers, padding=PADDING_ANG):
    """coords: list of (res_seq, x, y, z). residue_numbers: set of int.
    Returns [cx, cy, cz, sx, sy, sz] rounded to 2 decimals."""
    pts = np.array([(x, y, z) for seq, x, y, z in coords if seq in residue_numbers], dtype=float)
    if pts.size == 0:
        raise RuntimeError("No matching residues found in receptor for binding-site union")
    mn = pts.min(axis=0)
    mx = pts.max(axis=0)
    center = (mn + mx) / 2.0
    size = (mx - mn) + 2.0 * padding
    return [round(float(v), 2) for v in (center[0], center[1], center[2], size[0], size[1], size[2])]


def process_protein(task_root, protein):
    af3_dir = os.path.join(task_root, protein, "fine_screening", "AF3")
    summary_csv = os.path.join(af3_dir, "summary.csv")
    archives_dir = os.path.join(af3_dir, "output")

    if not os.path.isfile(summary_csv):
        return None, f"missing summary.csv"

    archives = sorted(
        f for f in os.listdir(archives_dir)
        if re.match(r"batch_\d+\.tar\.gz$", f)
    )
    if not archives:
        return None, "no batch archives"

    pdb_dir = os.path.join(task_root, "Input", "protein_file", protein)
    os.makedirs(pdb_dir, exist_ok=True)
    pdb_path = os.path.join(pdb_dir, f"{protein}.pdb")

    scratch = tempfile.mkdtemp(prefix=f"prep_{protein}_")
    try:
        cif_path = find_first_cif_in_archive(os.path.join(archives_dir, archives[0]), scratch)
        if cif_path is None:
            return None, f"no _model.cif in {archives[0]}"
        coords = cif_chain_a_to_pdb(cif_path, pdb_path)
    finally:
        shutil.rmtree(scratch, ignore_errors=True)

    res_nums = union_residue_numbers(summary_csv)
    if not res_nums:
        return None, "no binding_residues in summary.csv"

    box = compute_box(coords, res_nums)
    return box, None


def update_sequences_csv(task_root, boxes):
    seq_csv = os.path.join(task_root, "Input", "sequences.csv")
    df = pd.read_csv(seq_csv)
    if "PDB_ID" not in df.columns:
        df["PDB_ID"] = ""
    if "docking box" not in df.columns:
        df["docking box"] = ""

    name_col = "name" if "name" in df.columns else "protein_name"
    for protein, box in boxes.items():
        mask = df[name_col] == protein
        if not mask.any():
            print(f"  WARN: {protein} not in sequences.csv", flush=True)
            continue
        df.loc[mask, "PDB_ID"] = ""
        df.loc[mask, "docking box"] = json.dumps([box])

    df.to_csv(seq_csv, index=False)
    return seq_csv


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--task-root", default="/shared/cuteness_woundcare")
    p.add_argument("--only", nargs="*", default=None,
                   help="Restrict to these protein names")
    args = p.parse_args()

    proteins = discover_proteins(args.task_root)
    if args.only:
        proteins = [p for p in proteins if p in set(args.only)]

    print(f"Processing {len(proteins)} proteins", flush=True)
    boxes = {}
    failures = []
    for protein in proteins:
        try:
            box, err = process_protein(args.task_root, protein)
        except Exception as e:
            failures.append((protein, repr(e)))
            print(f"  {protein}: ERROR {e}", flush=True)
            continue
        if box is None:
            failures.append((protein, err))
            print(f"  {protein}: SKIP {err}", flush=True)
            continue
        boxes[protein] = box
        print(f"  {protein}: box {box}", flush=True)

    if boxes:
        seq_csv = update_sequences_csv(args.task_root, boxes)
        print(f"Updated {seq_csv} with {len(boxes)} entries", flush=True)
    else:
        print("No boxes produced; sequences.csv not modified", flush=True)

    if failures:
        print(f"\n{len(failures)} failures:", flush=True)
        for name, msg in failures:
            print(f"  {name}: {msg}", flush=True)


if __name__ == "__main__":
    main()
