###################################
# Write input files for AF3 with compounds
###################################

import pandas as pd 
import numpy as np
import json
import os
import glob
import argparse

def update_json(target_json, compound_smiles, output_path, i):
    # Step 1: Read the JSON file
    with open(target_json, "r") as f:
        data = json.load(f)  # Load JSON data as a Python dictionary

    data['name'] = data['name'] + str(i)
    data['version'] = 2
    # Step 2: Add a new argument (key-value pair)
    data["sequences"].append({})
    data["sequences"][1]["ligand"] = {}
    data["sequences"][1]["ligand"]["id"] = "Z"
    data["sequences"][1]["ligand"]["smiles"] = compound_smiles

    # Step 3: Save the updated JSON file
    with open(output_path, "w") as f:
        json.dump(data, f, indent=4)  # Save with indentation for readability
    print(output_path)
    return None

"""  -------- code used to copy file from previous AF3 output to this folder -----------
import shutil
data_json_list = glob.glob('/home/pacayang/Documents/B3/stage1_ML/preliminary/af3/control_input/**/*_data.json', recursive=True)
for file in data_json_list:
    outname = file.split("/")[-1]
    outpath = os.path.join("/home/pacayang/Documents/B3/stage2_AF3/input/protein_json_input", outname)
    shutil.copy2(file, outpath) 
"""
def main():
    parser = argparse.ArgumentParser(description="program to create AF3 input json", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--output-dir', help="destination for the output json files")
    parser.add_argument('--input-json', help="the pre-folded json file with MSA and templates")
    parser.add_argument('--smiles-file', help="csv file containing the smiles")
    parser.add_argument('--smiles-col', default='SMILES', help="the name of column containing the smiles")

    args = parser.parse_args()

    output_folder_gradparent = args.output_dir 
    smiles_list = pd.read_csv(args.smiles_file)[args.smiles_col]
    file = args.input_json
    target = (file.split('/')[-1]).split("_data.json")[0]
    i = 0 # track number of file written
    for smiles in smiles_list:
        filename = target + '_' + str(i) + ".json"
        output = os.path.join(args.output_dir, filename)
        update_json(file, smiles, output, i)
        i += 1
        if i % 1000 == 0:
            print("*")

if __name__ == "__main__":
    main()
