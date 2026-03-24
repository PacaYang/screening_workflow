import numpy as np
import pandas as pd
import os
import glob
import re
import argparse
from tqdm import tqdm

def collect_results(path):
    # Recursive glob to find all docking_affinities.txt files
    affinity_file_list = glob.glob(f'{path}/**/docking_affinities.txt', recursive=True)
    folder_lst = []
    idx_lst = []
    box_lst = []
    affinity_lst = []

    for file in affinity_file_list:
        # Parse path components relative to the base path
        rel = os.path.relpath(file, path)
        parts = rel.split(os.sep)

        # Expected layouts:
        #   new: {part}/lig{i}/box{j}/docking_affinities.txt  (4 parts)
        #   old: {part}/lig{i}/docking_affinities.txt          (3 parts)
        if len(parts) >= 4:
            folder = parts[-4].split('_')[-1]
            idx = parts[-3].replace('lig', '')
            box_match = re.match(r'box(\d+)', parts[-2])
            box = box_match.group(1) if box_match else '0'
        elif len(parts) >= 3:
            folder = parts[-3].split('_')[-1]
            idx = parts[-2].replace('lig', '')
            box = '0'
        else:
            # Fallback: try original parsing
            folder = (file.split('/')[-3]).split('_')[-1]
            idx = (file.split('/')[-2]).split('lig')[-1]
            box = '0'

        folder_lst.append(folder)
        idx_lst.append(idx)
        box_lst.append(box)
        with open(file, "r") as f:
            affinity = f.readline().strip()
        affinity_lst.append(affinity)

    tmp_df = pd.DataFrame(data={'folder':folder_lst, 'idx':idx_lst, 'box':box_lst, 'affinity':affinity_lst})
    return tmp_df

def process_results(scores_df, input_folder):
    df_dicts = {}

    smilesfile_list = glob.glob(f"{input_folder}/*.csv")
    for i in range(0,len(smilesfile_list)):
        df_dicts[str(i)] = pd.read_csv(f"{input_folder}/input_part_{i}.csv")
    
    def find_SMILES(input):
        try:
            smiles = df_dicts[str(input[0])].iloc[int(input[1])]['ligand_description']
        except:
            smiles = df_dicts[str(input[0])].iloc[int(input[1])]['SMILES']
        return smiles

    scores_df['SMILES'] = scores_df[['folder', 'idx']].apply(lambda x: find_SMILES(x), axis=1)
    return scores_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Summarize the scores for Vina prediction", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--vina-results-folder", type=str, help="Path of the folder containing all the Vina predicted results")
    parser.add_argument("--output-dir", type=str, help="Path to save the summary")
    parser.add_argument("--input-dir", type=str, help="Path to the input files. Used to map the results")

    args = parser.parse_args()
    tmp_df = collect_results(args.vina_results_folder)
    df = process_results(tmp_df, args.input_dir)
    df.to_csv(os.path.join(args.output_dir, 'results.csv'))
