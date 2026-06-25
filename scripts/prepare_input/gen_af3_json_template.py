###################################
# Write single-stage AF3 input JSONs using templates and no MSA
###################################
"""
Generate AlphaFold3 input JSONs for the template-based, single-stage flow.

Unlike the two-stage prefold path (gen_af3_json_protein.py ->
gen_af3_json_with_cmpds.py), this skips the data pipeline entirely:

  * Each protein chain gets ``unpairedMsa: ""`` and ``pairedMsa: ""`` (MSA-free).
  * Each protein chain gets a single-chain mmCIF template via the ``templates``
    list (queryIndices/templateIndices default to a full 1:1 residue mapping).
  * Multiple protein chains come from the comma-separated ``sequence`` cell.
  * N compounds are co-folded per prediction; the screened compound is id ``Z``.

The batch runner (run_af3_batch.sh) must be invoked with --norun_data_pipeline
(it already is) so AF3 consumes these MSA-free inputs directly.

Template file convention (single chain each, required by AF3):
    <template-dir>/<name>_<CHAIN>.cif    e.g. MCP_A.cif, MCP_B.cif
"""

import os
import json
import argparse

import pandas as pd

from seq_utils import parse_chains, chain_ids, ligand_ids


def build_protein_entry(chain_id, sequence, template_cif):
    """One AF3 protein sequence entry: MSA-free, optionally with a template."""
    protein = {
        "id": chain_id,
        "sequence": sequence,
        # Both must be set together to run MSA-free.
        "unpairedMsa": "",
        "pairedMsa": "",
        "templates": [],
    }
    if template_cif and os.path.isfile(template_cif):
        # Full 1:1 residue mapping; assumes the CIF chain matches the sequence.
        idx = list(range(len(sequence)))
        protein["templates"] = [{
            "mmcif": template_cif,
            "queryIndices": idx,
            "templateIndices": idx,
        }]
    return {"protein": protein}


def write_json(name, template_name, chain_seqs, ligand_smiles, template_dir, output_path, chunk_idx):
    cids = chain_ids(len(chain_seqs))
    data = {
        "name": f"{name}_{chunk_idx}",
        "dialect": "alphafold3",
        "version": 2,
        "modelSeeds": [99],
        "sequences": [],
    }

    for cid, seq in zip(cids, chain_seqs):
        template_cif = None
        if template_dir:
            # Template files follow the original-case convention: <name>_<CHAIN>.cif.
            template_cif = os.path.join(template_dir, f"{template_name}_{cid}.cif")
        data["sequences"].append(build_protein_entry(cid, seq, template_cif))

    for lid, smi in zip(ligand_ids(len(ligand_smiles)), ligand_smiles):
        data["sequences"].append({"ligand": {"id": lid, "smiles": smi}})

    with open(output_path, "w") as f:
        json.dump(data, f, indent=4)
    print(output_path)


def main():
    parser = argparse.ArgumentParser(
        description="Create single-stage AF3 input JSONs (template + no MSA, multi-chain/multi-ligand)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('--output-dir', required=True, help="destination for the output JSON files")
    parser.add_argument('--input-csv', required=True, help="sequences.csv with name + sequence columns")
    parser.add_argument('--name-col', default='name', help="column holding the protein/target name")
    parser.add_argument('--seq-col', default='sequence', help="column holding the (comma-separated) sequence(s)")
    parser.add_argument('--protein-name', required=True, help="target name to write inputs for")
    parser.add_argument('--smiles-file', required=True, help="CSV file containing the SMILES")
    parser.add_argument('--smiles-col', default='SMILES', help="column name for the SMILES")
    parser.add_argument('--template-dir', default=None,
                        help="dir with per-chain single-chain mmCIFs named <name>_<CHAIN>.cif")
    parser.add_argument('--n-ligands', type=int, default=1,
                        help="compounds to co-fold per prediction (default 1); screened compound is Z")

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    seq_df = pd.read_csv(args.input_csv)
    if args.protein_name not in seq_df[args.name_col].tolist():
        raise NameError(f"No protein named {args.protein_name} in {args.input_csv}")
    cell = seq_df.loc[seq_df[args.name_col] == args.protein_name, args.seq_col].iloc[0]
    chain_seqs = parse_chains(cell)

    smiles_list = pd.read_csv(args.smiles_file)[args.smiles_col].tolist()
    n = max(1, args.n_ligands)
    name_lower = args.protein_name.lower()

    for chunk_idx, start in enumerate(range(0, len(smiles_list), n)):
        ligands = smiles_list[start:start + n]
        output = os.path.join(args.output_dir, f"{name_lower}_{chunk_idx}.json")
        write_json(name_lower, args.protein_name, chain_seqs, ligands, args.template_dir, output, chunk_idx)
        if (chunk_idx + 1) % 1000 == 0:
            print("*")


if __name__ == "__main__":
    main()
