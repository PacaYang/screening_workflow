import json
import os
import pandas as pd
import numpy as np
import argparse
from tqdm import tqdm
import sys

sys.path.insert(0, os.path.dirname(__file__))
from binding_site_utils import compute_binding_site_cif

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

def gen_scores_df(path, include_binding_sites=False, cutoff_ang=10.0):
    # get target directory
    # loop over path
    columns = ['SMILES', 'chain_iptm', 'chain_pair_iptm', 'chain_pair_pae_min', 'chain_ptm', 'iptm', 'ptm', 'ranking_score', 'plddt', 'folder']
    if include_binding_sites:
        columns.append('binding_residues')

    # Collect compound dirs: handle flat compoundN/ and nested job_N/compoundN/ layouts.
    # A directory is a compound dir iff it contains <dirname>_data.json.
    def _is_compound_dir(name, full):
        return os.path.isdir(full) and 'token' not in name and \
            os.path.isfile(os.path.join(full, name + '_data.json'))

    candidates = []
    for entry in os.listdir(path):
        entry_path = os.path.join(path, entry)
        if not os.path.isdir(entry_path) or 'token' in entry:
            continue
        if _is_compound_dir(entry, entry_path):
            candidates.append((entry, entry_path))
        else:
            # Nested layout: job_N/compoundN/
            for sub in os.listdir(entry_path):
                sub_path = os.path.join(entry_path, sub)
                if _is_compound_dir(sub, sub_path):
                    candidates.append((sub, sub_path))

    df = None
    for folder, first_level_path in tqdm(candidates):
        dataFile = os.path.join(first_level_path, folder + '_data.json')
        with open(dataFile, "r") as file:
            data = json.load(file)
        smiles = data['sequences'][-1]['ligand']['smiles']
        summaryFile = os.path.join(first_level_path, folder + '_summary_confidences.json')
        plddtFile = os.path.join(first_level_path, folder + '_confidences.json')

        if not os.path.exists(summaryFile) or not os.path.exists(plddtFile):
            print("skip folder {}".format(first_level_path))
            continue

        scores = read_summary_confidences(summaryFile, plddtFile)
        scores.insert(0, smiles)
        scores.insert(9, folder)
        if include_binding_sites:
            cif_file = os.path.join(first_level_path, folder + '_model.cif')
            binding_residues = compute_binding_site_cif(
                cif_file, ligand_chain='Z', cutoff_ang=cutoff_ang
            )
            scores.append(binding_residues)
        if df is None:
            df = pd.DataFrame(data=[scores], columns=columns)
        else:
            df = pd.concat([df, pd.DataFrame([scores], columns=columns)], ignore_index=True)

    if df is None:
        df = pd.DataFrame(columns=columns)
    return df

def analyze(input_folder, output_dir, include_binding_sites=False, cutoff_ang=10.0):
    df = gen_scores_df(input_folder, include_binding_sites=include_binding_sites, cutoff_ang=cutoff_ang)
    output_name = os.path.join(output_dir, "summary.csv")
    df.to_csv(output_name)
    return None

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Summarize the scores for AF3 prediction", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--af3-results-folder", type=str, help="Path of the folder containing all the AF3 predicted results")
    parser.add_argument("--output-dir", type=str, help="Path to save the summary")
    parser.add_argument("--binding-sites", action="store_true", default=False,
                        help="Compute binding site residues from structure files (optional)")
    parser.add_argument("--distance-threshold", type=float, default=10.0,
                        help="Distance cutoff in Angstroms for binding site detection (default 10.0 = 1 nm)")
    args = parser.parse_args()

    analyze(args.af3_results_folder, args.output_dir,
            include_binding_sites=args.binding_sites,
            cutoff_ang=args.distance_threshold)

                          
