"""Backfill binding_site_center / binding_site_residues / num_binding_residues
into an existing RoseTTAFold summary.csv.

Usage:
    python enrich_rosettafold_summary.py --method-dir /path/to/fine_screening/RoseTTAFold
"""
import argparse
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from glob import glob

import pandas as pd
from tqdm import tqdm

import numpy as np
from Bio.PDB import PDBParser

def extract_binding_site_residues(pdb_path, distance_cutoff=5.0):
    """Copied from rosettafold_scores.py to avoid top-level torch import."""
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('complex', pdb_path)
        if 'A' not in structure[0] or 'B' not in structure[0]:
            return "", "", 0
        protein_chain = structure[0]['A']
        ligand_chain = structure[0]['B']
        ligand_residues = [res for res in ligand_chain if res.get_resname() == 'LG1']
        if not ligand_residues:
            return "", "", 0
        ligand_atoms = list(ligand_residues[0].get_atoms())
        ligand_coords = np.array([atom.coord for atom in ligand_atoms])
        ligand_center = ligand_coords.mean(axis=0)
        binding_residues = []
        for residue in protein_chain:
            if residue.id[0] != ' ':
                continue
            min_dist = min(
                np.linalg.norm(atom.coord - lig_atom.coord)
                for atom in residue.get_atoms()
                for lig_atom in ligand_atoms
            )
            if min_dist <= distance_cutoff:
                binding_residues.append(f"A:{residue.get_resname()}{residue.get_id()[1]}")
        center_str = f"{ligand_center[0]:.2f},{ligand_center[1]:.2f},{ligand_center[2]:.2f}"
        return ','.join(sorted(binding_residues)), center_str, len(binding_residues)
    except Exception as e:
        print(f"Warning: Could not extract binding site from {pdb_path}: {e}")
        return "", "", 0


def resolve_pdb(method_dir: str, folder: str) -> str | None:
    pattern = os.path.join(method_dir, "protein_ligand", "extracted", folder, "*.pdb")
    pdbs = glob(pattern)
    return pdbs[0] if pdbs else None


def process_row(args):
    method_dir, folder = args
    pdb = resolve_pdb(method_dir, folder)
    if pdb is None:
        return folder, "", "", 0
    residues, center, count = extract_binding_site_residues(pdb)
    return folder, residues, center, count


def main():
    parser = argparse.ArgumentParser(
        description="Enrich RoseTTAFold summary.csv with binding site coordinates",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--method-dir", required=True,
                        help="Path to fine_screening/RoseTTAFold/ directory")
    parser.add_argument("--workers", type=int, default=os.cpu_count(),
                        help="Parallel worker processes")
    parser.add_argument("--force", action="store_true",
                        help="Re-process rows that already have binding_site_center")
    args = parser.parse_args()

    summary_path = os.path.join(args.method_dir, "summary.csv")
    if not os.path.exists(summary_path):
        print(f"Error: {summary_path} not found")
        sys.exit(1)

    df = pd.read_csv(summary_path)

    for col in ("binding_site_residues", "binding_site_center", "num_binding_residues"):
        if col not in df.columns:
            df[col] = "" if col != "num_binding_residues" else 0

    if args.force:
        to_process = df.index.tolist()
    else:
        to_process = df.index[df["binding_site_center"].isna() | (df["binding_site_center"] == "")].tolist()

    print(f"Rows to process: {len(to_process)} / {len(df)}")
    if not to_process:
        print("Nothing to do.")
        return

    tasks = [(args.method_dir, str(df.at[i, "folder"])) for i in to_process]

    results = {}
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(process_row, t): t[1] for t in tasks}
        for fut in tqdm(as_completed(futures), total=len(futures), desc="Enriching RoseTTAFold"):
            folder, residues, center, count = fut.result()
            results[folder] = (residues, center, count)

    for i in to_process:
        folder = str(df.at[i, "folder"])
        if folder in results:
            residues, center, count = results[folder]
            df.at[i, "binding_site_residues"] = residues
            df.at[i, "binding_site_center"] = center
            df.at[i, "num_binding_residues"] = count

    tmp_path = summary_path + ".tmp"
    df.to_csv(tmp_path, index=False)
    os.replace(tmp_path, summary_path)
    enriched = df["binding_site_center"].notna() & (df["binding_site_center"] != "")
    print(f"Done. {enriched.sum()} / {len(df)} rows have binding_site_center.")


if __name__ == "__main__":
    main()
