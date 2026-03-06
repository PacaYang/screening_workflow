import json
import os
import pandas as pd
import numpy as np
from glob import glob
from pathlib import Path
from tqdm import tqdm 
import argparse

def get_affinity(affinity_file):
    with open(affinity_file, 'r') as f:
        data = json.load(f)
    f.close()

    entries = []
    for key, value in data.items():
        if key.startswith("affinity_pred_value"):
            suffix = key.replace("affinity_pred_value", "")
            prob_key = f"affinity_probability_binary{suffix}"
            probability = data.get(prob_key)
            entries.append((value, probability))

    # Find the entry with the lowest affinity_pred_value
    lowest = min(entries, key=lambda x: x[0])
    lowest_affinity, corresponding_probability = lowest
    return lowest_affinity, corresponding_probability

def read_summary_confidences(summaryFile, affinity_file):
    f = open(summaryFile, 'r')
    scores = json.load(f)
    f.close()
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

def gen_scores_df(path):
    # loop over path
    columns=['folder', 'SMILES', 'chain_iptm', 'pair_pde', 'chain_ptm', 'iptm', 'ptm', 'ranking_score', 'iplddt', 'affinity', 'probability']
    rows = []
    for folder in tqdm(os.listdir(path)):
        folder_path = os.path.join(path, folder)
        if not os.path.isdir(folder_path):
            continue

        # Read SMILES from smiles.txt
        smiles_file = os.path.join(folder_path, 'smiles.txt')
        smiles = ""
        if os.path.exists(smiles_file):
            with open(smiles_file, 'r') as f:
                smiles = f.read().strip()

        affinity_file = find_file(folder_path, 'affinity_*.json')
        summary_file = find_file(folder_path, 'confidence*.json')

        if affinity_file and summary_file:
            scores = read_summary_confidences(summary_file, affinity_file)
            scores.insert(0, smiles)
            scores.insert(0, folder)
            rows.append(scores)
        else:
            print(f"📁 Depth 1: {folder} is a file")

    df = pd.DataFrame(data=rows, columns=columns)
    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Summarize the scores for Boltz prediction", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--boltz-results-folder", type=str, help="Path of the folder containing all the Boltz predicted results")
    parser.add_argument("--output-dir", type=str, help="Path to save the summary")

    args = parser.parse_args()

    df = gen_scores_df(args.boltz_results_folder)
    output_name = os.path.join(args.output_dir, "summary.csv")
    df.to_csv(output_name, index=False)
