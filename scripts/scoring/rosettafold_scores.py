import os
import re
import io
import tarfile
import argparse
import torch
import yaml
import pandas as pd
from glob import glob
from tqdm import tqdm


def process_extracted(extracted_dir, protein_name):
    """Process scores from an extracted directory on disk."""
    rows = []

    # Find all _aux.pt files
    aux_files = glob(os.path.join(extracted_dir, "compound_*", "*_aux.pt"))
    print(f"Found {len(aux_files)} aux.pt files")

    for aux_path in tqdm(aux_files, desc="Processing compounds"):
        match = re.search(r'compound_(\d+)/', aux_path)
        if not match:
            continue
        cid = match.group(1)
        folder = f"compound_{cid}"
        compound_dir = os.path.dirname(aux_path)

        try:
            data = torch.load(aux_path, map_location='cpu', weights_only=False)
        except Exception as e:
            print(f"Warning: could not load {aux_path}: {e}")
            continue

        # Read SMILES from smiles.txt
        smiles = ""
        smiles_file = os.path.join(compound_dir, 'smiles.txt')
        if os.path.exists(smiles_file):
            with open(smiles_file, 'r') as f:
                smiles = f.read().strip()

        rows.append({
            'folder': folder,
            'SMILES': smiles,
            'mean_plddt': data.get('mean_plddt', None),
            'mean_pae': data.get('mean_pae', None),
            'pae_prot': data.get('pae_prot', None),
            'pae_inter': data.get('pae_inter', None),
        })

    return rows


def process_tar(tar_path, protein_name):
    """Extract scores and SMILES from a single RoseTTAFold tar.gz archive."""
    rows = []
    with tarfile.open(tar_path, 'r:gz') as tar:
        members = tar.getmembers()

        yaml_map = {}
        for m in members:
            match = re.search(r'config_compound_(\d+)\.yaml$', m.name)
            if match:
                yaml_map[match.group(1)] = m

        for m in members:
            match = re.search(r'compound_(\d+)/.*_aux\.pt$', m.name)
            if not match:
                continue
            cid = match.group(1)
            folder = f"compound_{cid}"

            try:
                f = tar.extractfile(m)
                data = torch.load(io.BytesIO(f.read()), map_location='cpu', weights_only=False)
            except Exception as e:
                print(f"Warning: could not load {m.name}: {e}")
                continue

            smiles = ""
            if cid in yaml_map:
                try:
                    yf = tar.extractfile(yaml_map[cid])
                    cfg = yaml.safe_load(yf)
                    smiles = cfg.get('sm_inputs', {}).get('B', {}).get('input', '')
                except Exception:
                    pass

            rows.append({
                'folder': folder,
                'SMILES': smiles,
                'mean_plddt': data.get('mean_plddt', None),
                'mean_pae': data.get('mean_pae', None),
                'pae_prot': data.get('pae_prot', None),
                'pae_inter': data.get('pae_inter', None),
            })

    return rows


def analyze(results_folder, protein_name, output_dir):
    # Check if there's an extracted directory with compound folders
    compound_dirs = glob(os.path.join(results_folder, "compound_*"))
    if compound_dirs:
        print(f"Found extracted compound directories in {results_folder}")
        all_rows = process_extracted(results_folder, protein_name)
    else:
        tar_files = sorted([
            os.path.join(results_folder, f)
            for f in os.listdir(results_folder)
            if f.endswith('.tar.gz')
        ])
        print(f"Found {len(tar_files)} tar.gz archives")
        all_rows = []
        for tar_path in tqdm(tar_files, desc="Processing archives"):
            all_rows.extend(process_tar(tar_path, protein_name))

    df = pd.DataFrame(all_rows)
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, "summary.csv")
    df.to_csv(output_path, index=False)
    print(f"Wrote {len(df)} rows to {output_path}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Summarize scores for RoseTTAFold predictions",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--rfaa-results-folder", type=str, required=True,
                        help="Path to folder containing RoseTTAFold results (extracted or tar.gz)")
    parser.add_argument("--protein-name", type=str, required=True,
                        help="Protein target name (e.g. IL4RA)")
    parser.add_argument("--output-dir", type=str, required=True,
                        help="Path to save the summary CSV")
    args = parser.parse_args()

    analyze(args.rfaa_results_folder, args.protein_name, args.output_dir)
