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

from seq_utils import parse_chains, chain_ids, ligand_ids


def write_boltz_template_input(output, chain_seqs, ligand_smiles, template_cif):
    """Write a single-stage Boltz2 YAML using a template and no MSA.

    Args:
        output: destination YAML path.
        chain_seqs: list of protein chain sequences (chains A, B, C, ...).
        ligand_smiles: list of ligand SMILES. The first is the screened
            compound and is assigned chain id ``Z`` so scorers keep working.
        template_cif: path to a (multi-chain) mmCIF template, or None.
    """
    cids = chain_ids(len(chain_seqs))
    sequences = []
    for cid, seq in zip(cids, chain_seqs):
        # ``msa: empty`` forces single-sequence mode (skip MSA search).
        sequences.append({"protein": {"id": cid, "sequence": seq, "msa": "empty"}})

    lids = ligand_ids(len(ligand_smiles))
    for lid, smi in zip(lids, ligand_smiles):
        sequences.append({"ligand": {"id": lid, "smiles": smi}})

    boltz_input_dict = {
        "version": 1,
        "sequences": sequences,
        # Affinity is scored for the screened compound (always id Z).
        "properties": [{"affinity": {"binder": "Z"}}],
    }

    if template_cif:
        boltz_input_dict["templates"] = [
            {"cif": os.path.abspath(template_cif), "chain_id": cids}
        ]

    with open(output, "w") as file:
        yaml.dump(boltz_input_dict, file, default_flow_style=False)
    return None


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

    parser.add_argument('--template-mode', action="store_true", default=False,
                    help='single-stage mode: skip MSA (msa: empty) and use a structural template. '
                         'Reads chains from the comma-separated sequence cell of --protein-file.')

    parser.add_argument('--template-cif',
                    type=str, default=None,
                    help='path to the mmCIF template (multi-chain OK). Used with --template-mode.')

    parser.add_argument('--n-ligands',
                    type=int, default=1,
                    help='number of compounds to co-fold per prediction (default 1). '
                         'selected.csv is chunked into groups of this size; the screened '
                         'compound is always ligand Z.')

    args = parser.parse_args()


    if args.template_mode:
        os.makedirs(args.output, exist_ok=True)
        protein_df = pd.read_csv(args.protein_file)
        cell = protein_df.loc[protein_df[args.name_col] == args.protein_name, args.seq_col].iloc[0]
        chain_seqs = parse_chains(cell)

        smiles_df = pd.read_csv(args.smiles_path)
        smiles_list = smiles_df['SMILES'].tolist()
        n = max(1, args.n_ligands)

        for chunk_idx, start in enumerate(range(0, len(smiles_list), n)):
            ligands = smiles_list[start:start + n]
            output_path = os.path.join(args.output, f"{chunk_idx}.yaml")
            write_boltz_template_input(output_path, chain_seqs, ligands, args.template_cif)

    elif args.protein_only:
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
