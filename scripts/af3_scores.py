import json
import os
import pandas as pd
import numpy as np
import argparse
from tqdm import tqdm 

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

def gen_scores_df(path):
    # get target directory
    # loop over path
    columns=['SMILES', 'chain_iptm', 'chain_pair_iptm', 'chain_pair_pae_min', 'chain_ptm', 'iptm', 'ptm', 'ranking_score', 'plddt','folder']
    for folder in tqdm(os.listdir(path)):
        first_level_path = os.path.join(path, folder)
        if os.path.isdir(first_level_path) and '_' not in folder and 'token' not in folder:  # Ensure it's a directory and skip any folder containing '_' in its name
            # print(f"📁 Depth 1: {first_level_path}")
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
            if 'df' not in locals():
                df = pd.DataFrame(data=[scores], columns=columns)
            else:
                df = pd.concat([df, pd.DataFrame([scores], columns=df.columns)], ignore_index=True)
        else:
            print(f"📁 Depth 1: {first_level_path} is a file")

    return df

def analyze(input_folder, output_dir):
    df = gen_scores_df(input_folder)
    output_name = os.path.join(output_dir, "summary.csv")
    df.to_csv(output_name)
    return None

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Summarize the scores for AF3 prediction", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--af3-results-folder", type=str, help="Path of the folder containing all the AF3 predicted results")
    parser.add_argument("--output-dir", type=str, help="Path to save the summary")
    args = parser.parse_args()

    analyze(args.af3_results_folder, args.output_dir)

                          
