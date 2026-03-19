import json
import os
import pandas as pd
import numpy as np
from pathlib import Path
from tqdm import tqdm
import argparse
from Bio.PDB import MMCIFParser

def extract_binding_site_residues(cif_path, distance_cutoff=5.0):
    """Extract binding site residues from Boltz2 CIF file.

    Returns:
        binding_site_residues: str - Comma-separated list (e.g., "A:CYS45,A:ARG67")
        binding_site_center: tuple - (x, y, z) coordinates
        num_binding_residues: int - Count of residues within cutoff
    """
    try:
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure('complex', cif_path)

        # Extract protein (Chain A) and ligand (Chain Z)
        if 'A' not in structure[0] or 'Z' not in structure[0]:
            return "", "", 0

        protein_chain = structure[0]['A']
        ligand_chain = structure[0]['Z']

        # Get ligand atoms (residue LIG1)
        ligand_atoms = [atom for residue in ligand_chain for atom in residue.get_atoms()]
        if not ligand_atoms:
            return "", "", 0

        ligand_coords = np.array([atom.coord for atom in ligand_atoms])
        ligand_center = ligand_coords.mean(axis=0)

        # Find binding site residues
        binding_residues = []
        for residue in protein_chain:
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
                binding_residues.append(f"A:{resname}{resid}")

        binding_site_str = ','.join(sorted(binding_residues))
        center_str = f"{ligand_center[0]:.2f},{ligand_center[1]:.2f},{ligand_center[2]:.2f}"
        return binding_site_str, center_str, len(binding_residues)
    except Exception as e:
        print(f"Warning: Could not extract binding site from {cif_path}: {e}")
        return "", "", 0

def get_affinity_entries(affinity_file):
    """Return all affinity/probability entries from a Boltz2 affinity JSON."""
    with open(affinity_file, 'r') as f:
        data = json.load(f)

    entries = []
    for key, value in data.items():
        if key.startswith("affinity_pred_value"):
            suffix = key.replace("affinity_pred_value", "")
            prob_key = f"affinity_probability_binary{suffix}"
            probability = data.get(prob_key)
            model_id = suffix if suffix else "0"
            entries.append(
                {
                    "model_id": model_id,
                    "affinity": value,
                    "probability": probability,
                }
            )

    if not entries:
        raise ValueError(f"No affinity_pred_value keys found in {affinity_file}")

    return sorted(entries, key=lambda x: float(x["affinity"]))

def get_affinity(affinity_file):
    entries = get_affinity_entries(affinity_file)
    best_entry = entries[0]
    return best_entry["affinity"], best_entry["probability"]

def read_summary_confidences(summary_file, affinity_file):
    with open(summary_file, 'r') as f:
        scores = json.load(f)

    lowest_affinity, corresponding_probability = get_affinity(affinity_file)
    chain_iptm = scores['pair_chains_iptm']['1']['0']
    pair_pae = scores['complex_pde']
    plddt_score = scores['complex_iplddt']

    return [chain_iptm, pair_pae, scores['chains_ptm']['1'], scores['iptm'], \
            scores['ptm'], scores['confidence_score'], plddt_score, lowest_affinity, corresponding_probability]

def find_file(root_folder, pattern):
    for path in Path(root_folder).rglob(pattern):
        return str(path)  # Return first match
    return None

def load_collection_state(state_file):
    """Load collection state from JSON file."""
    if not os.path.exists(state_file):
        return {'collected_compounds': {}}
    with open(state_file, 'r') as f:
        return json.load(f)

def update_collection_state(state_file, compound_id, has_data, warning=None):
    """Update collection state for a compound."""
    state = load_collection_state(state_file)
    state['collected_compounds'][compound_id] = {
        'timestamp': pd.Timestamp.now(tz='UTC').isoformat(),
        'status': 'collected',
        'has_data': has_data
    }
    if warning:
        state['collected_compounds'][compound_id]['warning'] = warning
    state['last_collection'] = pd.Timestamp.now(tz='UTC').isoformat()

    with open(state_file, 'w') as f:
        json.dump(state, f, indent=2)

def get_boltz2_row(folder, smiles, affinity_file, summary_file=None, cif_file=None):
    """Build a summary row from the files actually present in a Boltz2 output folder."""
    affinity_entries = get_affinity_entries(affinity_file)
    best_entry = affinity_entries[0]

    # Populate confidence fields only when a confidence JSON exists.
    if summary_file and os.path.exists(summary_file):
        score_values = read_summary_confidences(summary_file, affinity_file)
    else:
        score_values = [
            np.nan,  # chain_iptm
            np.nan,  # pair_pde
            np.nan,  # chain_ptm
            np.nan,  # iptm
            np.nan,  # ptm
            np.nan,  # ranking_score
            np.nan,  # iplddt
            best_entry["affinity"],
            best_entry["probability"],
        ]

    row = {
        "folder": folder,
        "SMILES": smiles,
        "chain_iptm": score_values[0],
        "pair_pde": score_values[1],
        "chain_ptm": score_values[2],
        "iptm": score_values[3],
        "ptm": score_values[4],
        "ranking_score": score_values[5],
        "iplddt": score_values[6],
        "affinity": score_values[7],
        "probability": score_values[8],
        "best_model": best_entry["model_id"],
    }

    # Preserve per-model affinities/probabilities when they exist.
    for entry in affinity_entries:
        suffix = entry["model_id"]
        row[f"affinity_model_{suffix}"] = entry["affinity"]
        row[f"probability_model_{suffix}"] = entry["probability"]

    if cif_file and os.path.exists(cif_file):
        binding_site_residues, binding_site_center, num_binding_residues = extract_binding_site_residues(cif_file)
    else:
        binding_site_residues, binding_site_center, num_binding_residues = "", "", 0

    row["binding_site_residues"] = binding_site_residues
    row["binding_site_center"] = binding_site_center
    row["num_binding_residues"] = num_binding_residues
    return row

def gen_scores_df(path, state_file=None, incremental=False):
    rows = []

    # Load state if incremental mode
    collected_ids = set()
    if incremental and state_file:
        state = load_collection_state(state_file)
        collected_ids = set(state.get('collected_compounds', {}).keys())
        print(f"Incremental mode: {len(collected_ids)} compounds already collected")

    all_folders = [
        f for f in os.listdir(path)
        if os.path.isdir(os.path.join(path, f)) and f.startswith("job_")
    ]

    # Filter to uncollected if incremental
    folders_to_process = [f for f in all_folders if f not in collected_ids] if incremental else all_folders

    if incremental and not folders_to_process:
        print("No new compounds to collect")
        return None

    print(f"Processing {len(folders_to_process)} compounds (total: {len(all_folders)})")

    for folder in tqdm(folders_to_process):
        folder_path = os.path.join(path, folder)
        if not os.path.isdir(folder_path):
            continue

        smiles_file = os.path.join(folder_path, 'smiles.txt')
        smiles = ""
        if os.path.exists(smiles_file):
            with open(smiles_file, 'r') as f:
                smiles = f.read().strip()

        try:
            affinity_file = find_file(folder_path, 'affinity_*.json')
            summary_file = find_file(folder_path, 'confidence*.json')
            cif_file = find_file(folder_path, '*_model_0.cif')

            if affinity_file:
                rows.append(
                    get_boltz2_row(
                        folder=folder,
                        smiles=smiles,
                        affinity_file=affinity_file,
                        summary_file=summary_file,
                        cif_file=cif_file,
                    )
                )

                # Update state if incremental
                if incremental and state_file:
                    update_collection_state(state_file, folder, True)
            else:
                warning = f"Missing affinity_*.json in {folder}"
                print(f"Warning: {warning}")
                if incremental and state_file:
                    update_collection_state(state_file, folder, False, warning)

        except Exception as e:
            warning = f"Failed to process {folder}: {str(e)}"
            print(f"Warning: {warning}")
            if incremental and state_file:
                update_collection_state(state_file, folder, False, warning)

    df = pd.DataFrame(data=rows)
    preferred_columns = [
        'folder', 'SMILES', 'chain_iptm', 'pair_pde', 'chain_ptm', 'iptm',
        'ptm', 'ranking_score', 'iplddt', 'affinity', 'probability',
        'best_model', 'binding_site_residues', 'binding_site_center',
        'num_binding_residues'
    ]
    extra_columns = [col for col in df.columns if col not in preferred_columns]
    df = df[[col for col in preferred_columns if col in df.columns] + sorted(extra_columns)]
    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Summarize the scores for Boltz prediction", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--boltz-results-folder", type=str, help="Path of the folder containing all the Boltz predicted results")
    parser.add_argument("--output-dir", type=str, help="Path to save the summary")
    parser.add_argument("--incremental", action="store_true", help="Enable incremental collection mode")
    parser.add_argument("--state-file", type=str, help="Path to collection state JSON file")
    parser.add_argument("--append", action="store_true", help="Append to existing CSV instead of overwriting")

    args = parser.parse_args()

    df = gen_scores_df(args.boltz_results_folder, args.state_file, args.incremental)

    if df is None or len(df) == 0:
        print("No new data to write")
    else:
        os.makedirs(args.output_dir, exist_ok=True)
        output_name = os.path.join(args.output_dir, "summary.csv")
        if args.append and os.path.exists(output_name):
            print(f"Appending {len(df)} rows to existing CSV")
            df.to_csv(output_name, mode='a', header=False, index=False)
        else:
            print(f"Writing {len(df)} rows to new CSV")
            df.to_csv(output_name, index=False)
