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
    columns=['folder', 'chain_iptm', 'pair_pde', 'chain_ptm', 'iptm', 'ptm', 'ranking_score', 'iplddt', 'affinity', 'probability']
    for folder in tqdm(os.listdir(path)):
        affinity_file = find_file(os.path.join(path, folder), 'affinity_*.json')
        summary_file = find_file(os.path.join(path, folder), 'confidence*.json')
        if affinity_file and summary_file:
            scores = read_summary_confidences(summary_file, affinity_file)
            scores.insert(0, folder)
            if 'df' not in locals():
                df = pd.DataFrame(data=[scores], columns=columns)
            else:
                df = pd.concat([df, pd.DataFrame([scores], columns=df.columns)], ignore_index=True)
        else:
            print(f"📁 Depth 1: {folder} is a file")

    return df

def map_smiles(df0, smiles_file, smiles_col):
    '''
    df0: the summary df without SMILES
    smiles_file: the path to the SMILES file
    '''
    df0["index"] = df0['folder'].str.split('_').str[-1].astype(int)
    smiles_df = pd.read_csv(smiles_file)[[smiles_col]]
    df = pd.merge(df0, smiles_df, left_on='index', right_index=True, how='left')
    return df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Summarize the scores for Boltz prediction", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--boltz-results-folder", type=str, help="Path of the folder containing all the Boltz predicted results")
    parser.add_argument("--output-dir", type=str, help="Path to save the summary")
    parser.add_argument("--smiles", type=str, help="Path to smiles file. Used to mapping the results")
    parser.add_argument("--smiles-col", default="SMILES", type=str, help="Col name for the smiles")

    args = parser.parse_args()

    df_0 = gen_scores_df(args.boltz_results_folder)
    df = map_smiles(df_0, args.smiles, args.smiles_col)
    output_name = os.path.join(args.output_dir, "summary.csv")
    df.to_csv(output_name, index=False)
