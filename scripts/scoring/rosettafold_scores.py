import os
import re
import io
import tarfile
import tempfile
import argparse
import torch
import yaml
import pandas as pd
from glob import glob
from tqdm import tqdm
import sys

sys.path.insert(0, os.path.dirname(__file__))
from binding_site_utils import compute_binding_site_pdb


def process_extracted(extracted_dir, protein_name, include_binding_sites=False, cutoff_ang=10.0):
    """Process scores from an extracted directory on disk."""
    rows = []

    # Build yaml map: supports config_{cid}.yaml and config_compound_{cid}.yaml
    yaml_map = {}
    for yf in glob(os.path.join(extracted_dir, "config_*.yaml")):
        match = re.search(r'config_(?:compound_)?(\d+)\.yaml$', yf)
        if match:
            yaml_map[match.group(1)] = yf

    # Find all _aux.pt files: supports compound_{cid}/*_aux.pt and {cid}/*_aux.pt
    aux_files = glob(os.path.join(extracted_dir, "*", "*_aux.pt"))
    print(f"Found {len(aux_files)} aux.pt files, {len(yaml_map)} config YAMLs")

    for aux_path in tqdm(aux_files, desc="Processing compounds"):
        match = re.search(r'/(?:compound_)?(\d+)/[^/]*_aux\.pt$', aux_path)
        if not match:
            continue
        cid = match.group(1)
        folder = f"compound_{cid}"

        try:
            data = torch.load(aux_path, map_location='cpu', weights_only=False)
        except Exception as e:
            print(f"Warning: could not load {aux_path}: {e}")
            continue

        smiles = ""
        if cid in yaml_map:
            try:
                with open(yaml_map[cid]) as f:
                    cfg = yaml.safe_load(f)
                smiles = cfg.get('sm_inputs', {}).get('Z', {}).get('input', '')
            except Exception:
                pass

        row = {
            'folder': folder,
            'SMILES': smiles,
            'mean_plddt': data.get('mean_plddt', None),
            'mean_pae': data.get('mean_pae', None),
            'pae_prot': data.get('pae_prot', None),
            'pae_inter': data.get('pae_inter', None),
        }

        if include_binding_sites:
            compound_dir = os.path.join(extracted_dir, folder)
            pdb_candidates = glob(os.path.join(compound_dir, '*_0.pdb'))
            if pdb_candidates:
                row['binding_residues'] = compute_binding_site_pdb(
                    pdb_candidates[0], ligand_chain='B', cutoff_ang=cutoff_ang
                )
            else:
                row['binding_residues'] = ''

        rows.append(row)

    return rows


def process_tar(tar_path, protein_name, include_binding_sites=False, cutoff_ang=10.0):
    """Extract scores and SMILES from a single RoseTTAFold tar.gz archive."""
    rows = []
    with tarfile.open(tar_path, 'r:gz') as tar:
        members = tar.getmembers()

        yaml_map = {}
        for m in members:
            # Support both config_compound_{cid}.yaml and config_{cid}.yaml
            match = re.search(r'config_(?:compound_)?(\d+)\.yaml$', m.name)
            if match:
                yaml_map[match.group(1)] = m

        # Build a map from compound id to *.pdb member for binding sites
        pdb_map = {}
        if include_binding_sites:
            for m in members:
                match = re.search(r'(?:compound_)?(\d+)/.*_0\.pdb$', m.name)
                if match:
                    pdb_map[match.group(1)] = m

        for m in members:
            # Support both compound_{cid}/..._aux.pt and {cid}/..._aux.pt
            match = re.search(r'/(?:compound_)?(\d+)/[^/]*_aux\.pt$', m.name)
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
                    smiles = cfg.get('sm_inputs', {}).get('Z', {}).get('input', '')
                except Exception:
                    pass

            row = {
                'folder': folder,
                'SMILES': smiles,
                'mean_plddt': data.get('mean_plddt', None),
                'mean_pae': data.get('mean_pae', None),
                'pae_prot': data.get('pae_prot', None),
                'pae_inter': data.get('pae_inter', None),
            }

            if include_binding_sites:
                binding_residues = ''
                if cid in pdb_map:
                    try:
                        pdb_bytes = tar.extractfile(pdb_map[cid]).read()
                        with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp:
                            tmp.write(pdb_bytes)
                            tmp_path = tmp.name
                        binding_residues = compute_binding_site_pdb(
                            tmp_path, ligand_chain='B', cutoff_ang=cutoff_ang
                        )
                        os.unlink(tmp_path)
                    except Exception:
                        binding_residues = ''
                row['binding_residues'] = binding_residues

            rows.append(row)

    return rows


def analyze(results_folder, protein_name, output_dir, include_binding_sites=False, cutoff_ang=10.0):
    # Check if there's an extracted directory with compound or numeric id folders
    compound_dirs = glob(os.path.join(results_folder, "compound_*"))
    numeric_dirs = [
        d for d in glob(os.path.join(results_folder, "*"))
        if os.path.isdir(d) and os.path.basename(d).isdigit()
    ]
    if compound_dirs or numeric_dirs:
        print(f"Found extracted directories in {results_folder} (compound: {len(compound_dirs)}, numeric: {len(numeric_dirs)})")
        all_rows = process_extracted(results_folder, protein_name,
                                     include_binding_sites=include_binding_sites,
                                     cutoff_ang=cutoff_ang)
    else:
        tar_files = sorted([
            os.path.join(results_folder, f)
            for f in os.listdir(results_folder)
            if f.endswith('.tar.gz')
        ])
        print(f"Found {len(tar_files)} tar.gz archives")
        all_rows = []
        for tar_path in tqdm(tar_files, desc="Processing archives"):
            all_rows.extend(process_tar(tar_path, protein_name,
                                        include_binding_sites=include_binding_sites,
                                        cutoff_ang=cutoff_ang))

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
    parser.add_argument("--binding-sites", action="store_true", default=False,
                        help="Compute binding site residues from structure files (optional)")
    parser.add_argument("--distance-threshold", type=float, default=10.0,
                        help="Distance cutoff in Angstroms for binding site detection (default 10.0 = 1 nm)")
    args = parser.parse_args()

    analyze(args.rfaa_results_folder, args.protein_name, args.output_dir,
            include_binding_sites=args.binding_sites,
            cutoff_ang=args.distance_threshold)
