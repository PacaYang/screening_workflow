###################################
# Write input files for AF3 with 2 competing ligands (Z and Y)
###################################

import pandas as pd
import json
import os
import argparse


def update_json(target_json, smiles_z, smiles_y, output_path, i):
    with open(target_json, "r") as f:
        data = json.load(f)

    data['name'] = data['name'] + str(i)
    data['version'] = 2

    # Ligand Z: query/candidate compound
    data["sequences"].append({
        "ligand": {
            "id": "Z",
            "smiles": smiles_z
        }
    })

    # Ligand Y: reference/competitor compound
    data["sequences"].append({
        "ligand": {
            "id": "Y",
            "smiles": smiles_y
        }
    })

    with open(output_path, "w") as f:
        json.dump(data, f, indent=4)
    print(output_path)


def main():
    parser = argparse.ArgumentParser(
        description="Generate AF3 input JSONs for competitive folding (1 protein + 2 ligands)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--output-dir', required=True,
                        help="Destination for the output JSON files")
    parser.add_argument('--input-json', required=True,
                        help="Prefold _data.json with MSA and templates")
    parser.add_argument('--smiles-file', required=True,
                        help="CSV file containing query SMILES (ligand Z)")
    parser.add_argument('--smiles-col', default='SMILES',
                        help="Column name for SMILES in smiles-file")
    parser.add_argument('--reference-smiles', required=True,
                        help="SMILES string for the reference/competitor ligand Y")

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    smiles_list = pd.read_csv(args.smiles_file)[args.smiles_col].tolist()
    target = os.path.basename(args.input_json).split("_data.json")[0]

    for i, smiles_z in enumerate(smiles_list):
        filename = f"{target}_{i}.json"
        output = os.path.join(args.output_dir, filename)
        update_json(args.input_json, smiles_z, args.reference_smiles, output, i)
        if (i + 1) % 1000 == 0:
            print("*")


if __name__ == "__main__":
    main()
