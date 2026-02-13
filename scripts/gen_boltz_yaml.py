###################################
# Write input files for Boltz2
###################################

import json
import numpy as np
import pandas as pd
import yaml
import os
import argparse
from tqdm import tqdm


def write_boltz_input(output, protein_sequence, msa_file_path=None, ligand=None, ligand_type="ligand"):
    # currently only support 1 protien + 1 ligand. The ligand can be any type supported by boltz
    # output: output file name
    # ligand: protein/RNA/DNA sequences or SMILES

    # check ligand type first
    if ligand_type not in ['ligand','protein','dna','rna']:
        print("unsupported ligand type, check boltz documents")
        return None

    # write input 
    if msa_file_path:
        protein_dict = {"protein":{"id":"A", "sequence": protein_sequence, "msa": msa_file_path}}
    else:
        protein_dict = {"protein":{"id":"A", "sequence": protein_sequence}}
    
    if ligand:
        ligand_dict = {ligand_type:{"id":"Z", "smiles": ligand}}
        boltz_input_dict = {"version": 1, "sequences":[protein_dict, ligand_dict], "properties":[{"affinity":{"binder":"Z"}}]}
    
    else:
        boltz_input_dict = {"version": 1, "sequences":[protein_dict],}

    with open(output, 'w') as file:
        yaml.dump(boltz_input_dict, file, default_flow_style=False)
    file.close()
    return None

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Write boltz input.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--output",
                    type=str,
                    help='destination of written files')

    parser.add_argument("--msa",
                    type=str,
                    help='destination of msa csv')
    
    parser.add_argument("--sequence",
                    type=str,
                    help='protein sequence. If it is missing, protein file and protein name will be used to extract the sequence from the file.')

    parser.add_argument("--smiles-path",
                    type=str,
                    help='destination of smiles csv')

    parser.add_argument("--protein-file",
                    default=None,
                    type=str,
                    help='destination of csv file containing the protein names and sequences.')   

    parser.add_argument('--name-col', 
                    type=str, default='name', help="the name of column containing the name. Use with --protein-file.")
    
    parser.add_argument('--seq-col', 
                    type=str, default='sequence', help="the name of column containing the sequence. Use with --protein-file.")
    
    parser.add_argument('--protein-name', 
                    help='protein to write. Use with --protein-file.')

    parser.add_argument('--protein-only', action="store_true", default=False,
                    help='write inputs for protein only. Used to prefold the protein. Ligands will be ignored with this option.')

    args = parser.parse_args()

    
    if args.protein_only:
        protein_df = pd.read_csv(args.protein_file)
        seq = protein_df.loc[protein_df[args.name_col]==args.protein_name, args.seq_col].iloc[0]
        output_path = os.path.join(args.output, f"{args.protein_name}.yaml")
        write_boltz_input(output_path, seq)

    else:
        smiles_df = pd.read_csv(args.smiles_path)
        
        for i, smi in tqdm(enumerate(smiles_df['SMILES'])):
            output_path = os.path.join(args.output, f"{i}.yaml")
            if args.sequence:
                write_boltz_input(output_path, args.sequence, args.msa, smi, ligand_type="ligand")
            else:
                protein_df = pd.read_csv(args.protein_file)
                seq = protein_df.loc[protein_df[args.name_col]==args.protein_name, args.seq_col].iloc[0]
                write_boltz_input(output_path, seq, args.msa, smi, ligand_type="ligand")
