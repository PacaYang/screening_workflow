import os
import re
import io
import tarfile
import argparse
import torch
import yaml
import pandas as pd
import numpy as np
from glob import glob
from tqdm import tqdm
from Bio.PDB import PDBParser


def extract_binding_site_residues(pdb_path, distance_cutoff=5.0):
    """Extract binding site residues from a RoseTTAFold PDB file.

    Multichain-aware: RFAA relabels output chains sequentially (A, B, C... by config
    order, proteins first then ligands), so the ligand chain is NOT fixed at 'B'. We
    locate the ligand as the residue named 'LG1' wherever it appears, and treat every
    other chain as protein. For multi-ligand co-folds the first LG1 found (the primary
    screened compound, id Z) defines the binding site.

    Returns:
        binding_site_residues: str - Comma-separated list (e.g., "A:CYS45,B:ARG67")
        binding_site_center: str - "x,y,z" of the ligand centroid
        num_binding_residues: int - Count of residues within cutoff
    """
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('complex', pdb_path)
        model = structure[0]

        # Find the ligand residue (resname LG1) across all chains.
        ligand_residue = None
        ligand_chain_id = None
        for chain in model:
            for res in chain:
                if res.get_resname() == 'LG1':
                    ligand_residue = res
                    ligand_chain_id = chain.id
                    break
            if ligand_residue is not None:
                break

        if ligand_residue is None:
            return "", "", 0

        ligand_atoms = list(ligand_residue.get_atoms())
        ligand_coords = np.array([atom.coord for atom in ligand_atoms])
        ligand_center = ligand_coords.mean(axis=0)

        # Every chain that is not the ligand chain is protein.
        binding_residues = []
        for chain in model:
            if chain.id == ligand_chain_id:
                continue
            for residue in chain:
                if residue.id[0] != ' ':  # Skip hetero residues
                    continue
                min_dist = float('inf')
                for atom in residue.get_atoms():
                    for lig_atom in ligand_atoms:
                        dist = np.linalg.norm(atom.coord - lig_atom.coord)
                        min_dist = min(min_dist, dist)

                if min_dist <= distance_cutoff:
                    resname = residue.get_resname()
                    resid = residue.get_id()[1]
                    binding_residues.append(f"{chain.id}:{resname}{resid}")

        binding_site_str = ','.join(sorted(binding_residues))
        center_str = f"{ligand_center[0]:.2f},{ligand_center[1]:.2f},{ligand_center[2]:.2f}"
        return binding_site_str, center_str, len(binding_residues)
    except Exception as e:
        print(f"Warning: Could not extract binding site from {pdb_path}: {e}")
        return "", "", 0


def process_extracted(extracted_dir, protein_name):
    """Process scores from an extracted directory on disk."""
    rows = []

    # Find all _aux.pt files
    aux_files = glob(os.path.join(extracted_dir, "compound_*", "*_aux.pt"))
    print(f"Found {len(aux_files)} aux.pt files")

    for aux_path in tqdm(aux_files, desc="Processing compounds"):
        match = re.search(r'compound_(\d+)/', aux_path)
        if not match:
            continue
        cid = match.group(1)
        folder = f"compound_{cid}"
        compound_dir = os.path.dirname(aux_path)

        try:
            data = torch.load(aux_path, map_location='cpu', weights_only=False)
        except Exception as e:
            print(f"Warning: could not load {aux_path}: {e}")
            continue

        # Read SMILES from smiles.txt
        smiles = ""
        smiles_file = os.path.join(compound_dir, 'smiles.txt')
        if os.path.exists(smiles_file):
            with open(smiles_file, 'r') as f:
                smiles = f.read().strip()

        # Extract binding site information from PDB file
        pdb_file = glob(os.path.join(compound_dir, "*.pdb"))
        if pdb_file:
            binding_site_residues, binding_site_center, num_binding_residues = extract_binding_site_residues(pdb_file[0])
        else:
            binding_site_residues, binding_site_center, num_binding_residues = "", "", 0

        rows.append({
            'folder': folder,
            'SMILES': smiles,
            'mean_plddt': data.get('mean_plddt', None),
            'mean_pae': data.get('mean_pae', None),
            'pae_prot': data.get('pae_prot', None),
            'pae_inter': data.get('pae_inter', None),
            'binding_site_residues': binding_site_residues,
            'binding_site_center': binding_site_center,
            'num_binding_residues': num_binding_residues,
        })

    return rows


def process_tar(tar_path, protein_name):
    """Extract scores and SMILES from a single RoseTTAFold tar.gz archive."""
    rows = []
    with tarfile.open(tar_path, 'r:gz') as tar:
        members = tar.getmembers()

        yaml_map = {}
        pdb_map = {}
        for m in members:
            match = re.search(r'config_compound_(\d+)\.yaml$', m.name)
            if match:
                yaml_map[match.group(1)] = m
            match_pdb = re.search(r'compound_(\d+)/.*\.pdb$', m.name)
            if match_pdb:
                pdb_map[match_pdb.group(1)] = m

        for m in members:
            match = re.search(r'compound_(\d+)/.*_aux\.pt$', m.name)
            if not match:
                continue
            cid = match.group(1)
            folder = f"compound_{cid}"

            try:
                f = tar.extractfile(m)
                data = torch.load(io.BytesIO(f.read()), map_location='cpu', weights_only=False)
            except Exception as e:
                print(f"Warning: could not load {m.name}: {e}")
                continue

            smiles = ""
            if cid in yaml_map:
                try:
                    yf = tar.extractfile(yaml_map[cid])
                    cfg = yaml.safe_load(yf)
                    sm = cfg.get('sm_inputs', {}) or {}
                    # Primary ligand is id Z; fall back to the first sm input.
                    entry = sm.get('Z') or (next(iter(sm.values())) if sm else {})
                    smiles = str((entry or {}).get('input', '')).strip()
                except Exception:
                    pass

            # Extract binding site information from PDB file in tar
            binding_site_residues, binding_site_center, num_binding_residues = "", "", 0
            if cid in pdb_map:
                try:
                    pdb_f = tar.extractfile(pdb_map[cid])
                    pdb_content = pdb_f.read()
                    # Write to temporary file for BioPython parsing
                    import tempfile
                    with tempfile.NamedTemporaryFile(mode='wb', suffix='.pdb', delete=False) as tmp:
                        tmp.write(pdb_content)
                        tmp_path = tmp.name
                    binding_site_residues, binding_site_center, num_binding_residues = extract_binding_site_residues(tmp_path)
                    os.unlink(tmp_path)
                except Exception as e:
                    print(f"Warning: could not extract binding site from {pdb_map[cid].name}: {e}")

            rows.append({
                'folder': folder,
                'SMILES': smiles,
                'mean_plddt': data.get('mean_plddt', None),
                'mean_pae': data.get('mean_pae', None),
                'pae_prot': data.get('pae_prot', None),
                'pae_inter': data.get('pae_inter', None),
                'binding_site_residues': binding_site_residues,
                'binding_site_center': binding_site_center,
                'num_binding_residues': num_binding_residues,
            })

    return rows


def analyze(results_folder, protein_name, output_dir):
    # Check if there's an extracted directory with compound folders
    compound_dirs = glob(os.path.join(results_folder, "compound_*"))
    if compound_dirs:
        print(f"Found extracted compound directories in {results_folder}")
        all_rows = process_extracted(results_folder, protein_name)
    else:
        tar_files = sorted([
            os.path.join(results_folder, f)
            for f in os.listdir(results_folder)
            if f.endswith('.tar.gz')
        ])
        print(f"Found {len(tar_files)} tar.gz archives")
        all_rows = []
        for tar_path in tqdm(tar_files, desc="Processing archives"):
            all_rows.extend(process_tar(tar_path, protein_name))

    df = pd.DataFrame(all_rows)
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, "summary.csv")
    df.to_csv(output_path, index=False)
    print(f"Wrote {len(df)} rows to {output_path}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Summarize scores for RoseTTAFold predictions",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--rfaa-results-folder", type=str, required=True,
                        help="Path to folder containing RoseTTAFold results (extracted or tar.gz)")
    parser.add_argument("--protein-name", type=str, required=True,
                        help="Protein target name (e.g. IL4RA)")
    parser.add_argument("--output-dir", type=str, required=True,
                        help="Path to save the summary CSV")
    args = parser.parse_args()

    analyze(args.rfaa_results_folder, args.protein_name, args.output_dir)
