###################################
# Write input files for AF3 without compounds (proteins only)
###################################

import pandas as pd 
import numpy as np
import json
import os
import glob
import argparse

def write_json(seq, protein_name, output_dir):
    # Step 1: Read the JSON file
    data = {}

    data['name'] = protein_name
    data['dialect'] = 'alphafold3'
    data['version'] = 1
    data['sequences'] = []
    data['modelSeeds'] = [99]
    # Step 2: Add a new argument (key-value pair)
    data["sequences"].append({})
    data["sequences"][0]["protein"] = {}
    data["sequences"][0]["protein"]["id"] = "A"
    data["sequences"][0]["protein"]["sequence"] = seq

    # Step 3: Save the updated JSON file
    with open(os.path.join(output_dir,f"{protein_name}.json"), "w") as f:
        json.dump(data, f, indent=4)  # Save with indentation for readability

    return None

def main():
    parser = argparse.ArgumentParser(description="program to create AF3 input json with proteins only", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--output-dir', help="destination for the output json files")
    parser.add_argument('--input-csv', help="input csv file containing the name and the sequence of the protein")
    parser.add_argument('--name-col', default='name', help="csv file containing the smiles")
    parser.add_argument('--seq-col', default='sequence', help="the name of column containing the smiles")
    parser.add_argument('--protein-name', help='protein to write')

    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir, exist_ok=True)

    if args.protein_name in pd.read_csv(args.input_csv)[args.name_col].to_list():
        df = pd.read_csv(args.input_csv)
        seq = df.loc[df[args.name_col] == args.protein_name, args.seq_col].iloc[0]
        write_json(seq, args.protein_name, args.output_dir)
    else:
        raise NameError(f"No protein named {args.protein_name} in the {args.input_csv} file")

if __name__ == "__main__":
    main()
