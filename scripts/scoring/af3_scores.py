import json
import os
import pandas as pd
import numpy as np
import argparse
from tqdm import tqdm
from Bio.PDB import MMCIFParser 

def extract_binding_site_residues(cif_path, distance_cutoff=5.0):
    """Extract binding site residues from AF3 CIF file.

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

        # Get ligand atoms
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

def get_plddts(json_summary):
    def find_indices(lst, value):
        return [i for i, x in enumerate(lst) if x == value]
        # Open and load the JSON file
    with open(json_summary, "r") as file:
        data = json.load(file)
    file.close()
    indices = find_indices(data['atom_chain_ids'], 'Z')
    plddts_lst = []
    for i in indices:
        plddts_lst.append(data['atom_plddts'][i])

    return plddts_lst

def read_summary_confidences(summaryFile, plddtFile):
    f = open(summaryFile, 'r')
    scores = json.load(f)
    f.close()
    ligand_chain_id = str(len(scores['chain_iptm']) + 1)

    plddt_score = np.average(get_plddts(plddtFile))
    if None in scores['chain_pair_iptm'][-1]:
        pair_iptm = -1
    else:
        pair_iptm = max(scores['chain_pair_iptm'][-1][:-1])


    if None in scores['chain_pair_pae_min'][-1]:
        pair_pae = -1
    else:
        pair_pae = min(scores['chain_pair_pae_min'][-1][:-1])

    if scores['chain_iptm'][-1] == None:
        chain_iptm = -1
    else:
        chain_iptm = scores['chain_iptm'][-1]

    # print(pair_iptm, scores['chain_ptm'][-1])

    return [chain_iptm, pair_iptm, pair_pae, scores['chain_ptm'][-1], scores['iptm'], \
            scores['ptm'], scores['ranking_score'], plddt_score]

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

def gen_scores_df(path, state_file=None, incremental=False):
    # get target directory
    # loop over path
    columns=['SMILES', 'chain_iptm', 'chain_pair_iptm', 'chain_pair_pae_min', 'chain_ptm', 'iptm', 'ptm', 'ranking_score', 'plddt','folder', 'binding_site_residues', 'binding_site_center', 'num_binding_residues']
    rows = []

    # Load state if incremental mode
    collected_ids = set()
    if incremental and state_file:
        state = load_collection_state(state_file)
        collected_ids = set(state.get('collected_compounds', {}).keys())
        print(f"Incremental mode: {len(collected_ids)} compounds already collected")

    all_folders = [f for f in os.listdir(path) if os.path.isdir(os.path.join(path, f))
                   and 'token' not in f]

    # Filter to uncollected if incremental
    folders_to_process = [f for f in all_folders if f not in collected_ids] if incremental else all_folders

    if incremental and not folders_to_process:
        print("No new compounds to collect")
        return None

    print(f"Processing {len(folders_to_process)} compounds (total: {len(all_folders)})")

    for folder in tqdm(folders_to_process):
        first_level_path = os.path.join(path, folder)
        if os.path.isdir(first_level_path) and 'token' not in folder:
            # Read SMILES from smiles.txt
            smiles_file = os.path.join(first_level_path, 'smiles.txt')
            if not os.path.exists(smiles_file):
                warning = f"No smiles.txt found in {first_level_path}"
                print(f"Warning: {warning}")
                if incremental and state_file:
                    update_collection_state(state_file, folder, False, warning)
                continue

            try:
                with open(smiles_file, 'r') as f:
                    smiles = f.read().strip()

                # Detect file prefix from summary_confidences.json files in folder
                import glob as _glob
                summary_matches = _glob.glob(os.path.join(first_level_path, '*_summary_confidences.json'))
                if summary_matches:
                    prefix = os.path.basename(summary_matches[0]).replace('_summary_confidences.json', '')
                else:
                    prefix = folder

                summaryFile = os.path.join(first_level_path, prefix + '_summary_confidences.json')
                plddtFile = os.path.join(first_level_path, prefix + '_confidences.json')
                cifFile = os.path.join(first_level_path, prefix + '_model.cif')

                if not os.path.exists(summaryFile) or not os.path.exists(plddtFile):
                    warning = f"Missing required JSON files in {first_level_path}"
                    print(f"Warning: {warning}")
                    if incremental and state_file:
                        update_collection_state(state_file, folder, False, warning)
                    continue

                scores = read_summary_confidences(summaryFile, plddtFile)
                scores.insert(0, smiles)
                scores.insert(9, folder)

                # Extract binding site information
                if os.path.exists(cifFile):
                    binding_site_residues, binding_site_center, num_binding_residues = extract_binding_site_residues(cifFile)
                else:
                    binding_site_residues, binding_site_center, num_binding_residues = "", "", 0

                scores.extend([binding_site_residues, binding_site_center, num_binding_residues])
                rows.append(scores)

                # Update state if incremental
                if incremental and state_file:
                    update_collection_state(state_file, folder, True)

            except Exception as e:
                warning = f"Failed to process {folder}: {str(e)}"
                print(f"Warning: {warning}")
                if incremental and state_file:
                    update_collection_state(state_file, folder, False, warning)
        else:
            print(f"📁 Depth 1: {first_level_path} is a file")

    df = pd.DataFrame(data=rows, columns=columns)
    return df

def analyze(input_folder, output_dir, state_file=None, incremental=False, append=False):
    df = gen_scores_df(input_folder, state_file, incremental)

    if df is None or len(df) == 0:
        print("No new data to write")
        return None

    output_name = os.path.join(output_dir, "summary.csv")

    # Append or write
    if append and os.path.exists(output_name):
        print(f"Appending {len(df)} rows to existing CSV")
        df.to_csv(output_name, mode='a', header=False, index=False)
    else:
        print(f"Writing {len(df)} rows to new CSV")
        df.to_csv(output_name, index=False)

    return None

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Summarize the scores for AF3 prediction", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--af3-results-folder", type=str, help="Path of the folder containing all the AF3 predicted results")
    parser.add_argument("--output-dir", type=str, help="Path to save the summary")
    parser.add_argument("--incremental", action="store_true", help="Enable incremental collection mode")
    parser.add_argument("--state-file", type=str, help="Path to collection state JSON file")
    parser.add_argument("--append", action="store_true", help="Append to existing CSV instead of overwriting")
    args = parser.parse_args()

    analyze(args.af3_results_folder, args.output_dir, args.state_file, args.incremental, args.append)

                          
