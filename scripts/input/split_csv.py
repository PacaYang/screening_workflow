###################################
# split the csv file containing the SMILES from initial screening
# Used as input files for Diffdock and autodock Vina
###################################

import pandas as pd 
import os
import argparse
from tqdm import tqdm

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="script to separate csv files for Diffdock & Vina parallel computations. " \
    "The generated files will be used as input for Vina and DiffDock, as well for mapping the results for Vina.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--protein-name", type=str, help="the protein name")
    parser.add_argument("--pdb-file", type=str, help="the protein pdb file")
    parser.add_argument("--smiles-file", type=str, help="the SMILES file")
    parser.add_argument("--chunk-size", type=int, help="number of entries in one file")
    parser.add_argument("--output-dir", type=str, help="the output directory")
    parser.add_argument('--smiles-col', type=str, default='SMILES', help="the name of column containing the smiles.")

    args = parser.parse_args()

    smiles_list = pd.read_csv(args.smiles_file)[args.smiles_col]
    if len(smiles_list)%args.chunk_size != 0:
        num_file = len(smiles_list)//args.chunk_size + 1
    else:
        num_file = len(smiles_list)//args.chunk_size

    for i in tqdm(range(num_file-1)):
        comp_names = []
        for j in range(args.chunk_size):
             comp_names.append(args.protein_name + "_" + str(i*args.chunk_size + j))
        tmp_df = pd.DataFrame(data={"complex_name":comp_names, "protein_path": [args.pdb_file]*args.chunk_size, 
                  "ligand_description":smiles_list[i*args.chunk_size: (i+1)*args.chunk_size],})
        tmp_df['protein_sequence']=''
        output_name = os.path.join(args.output_dir, (f"input_part_{i}.csv"))
        tmp_df.to_csv(output_name, index=False)
    
    # write the last input file
    comp_names = []
    length_to_write = len(smiles_list)-(num_file-1)*args.chunk_size
    for j in range(length_to_write):
        comp_names.append(args.protein_name + "_" + str((num_file-1)*args.chunk_size + j))
    tmp_df = pd.DataFrame(data={"complex_name":comp_names, "protein_path": [args.pdb_file]*length_to_write, 
                  "ligand_description":smiles_list[(num_file-1)*args.chunk_size: (num_file-1)*args.chunk_size+length_to_write],})
    tmp_df['protein_sequence']=''
    output_name = os.path.join(args.output_dir, (f"input_part_{num_file-1}.csv"))
    tmp_df.to_csv(output_name, index=False)
    print(f"Finished generating the csv files. Saved to {args.output_dir}")


