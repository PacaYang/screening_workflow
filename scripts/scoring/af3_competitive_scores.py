"""
Score parser for AF3 competitive folding results (1 protein + 2 ligands: Z and Y).

Chain layout in output:
  index 0 → protein chain A
  index 1 → ligand Z (query/candidate compound)
  index 2 → ligand Y (reference/competitor compound)

Metrics reported per ligand:
  chain_iptm, pair_iptm_vs_protein, pair_pae_min_vs_protein, chain_ptm, plddt

Additional:
  iptm, ptm, ranking_score (whole-complex)
  binding_site_Z, binding_site_Y — protein residues within --distance-threshold Å of each ligand
"""

import json
import os
import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm


# ============================================================================
# CIF parser
# ============================================================================

def parse_atom_site(cif_path):
    """
    Parse the _atom_site loop from an mmCIF file.

    Returns a list of dicts, one per atom, with keys matching the column
    names that follow _atom_site. in the loop header.
    Only columns present in the file are returned; caller should handle
    missing keys gracefully.
    """
    cols = []
    rows = []
    state = 'searching'  # searching | in_loop_start | in_atom_site_header | in_atom_site_data | in_other_loop

    with open(cif_path) as f:
        for line in f:
            stripped = line.strip()

            if not stripped or stripped.startswith('#'):
                continue

            if state == 'searching':
                if stripped == 'loop_':
                    state = 'in_loop_start'
                    cols = []

            elif state == 'in_loop_start':
                if stripped.startswith('_atom_site.'):
                    cols.append(stripped[len('_atom_site.'):])
                    state = 'in_atom_site_header'
                elif stripped.startswith('_'):
                    state = 'in_other_loop'
                # else: data row of some other block, go back to searching
                else:
                    state = 'searching'

            elif state == 'in_atom_site_header':
                if stripped.startswith('_atom_site.'):
                    cols.append(stripped[len('_atom_site.'):])
                elif stripped == 'loop_':
                    # Another loop started before any data — atom_site was empty
                    break
                elif stripped.startswith('data_'):
                    break
                else:
                    # First data row
                    state = 'in_atom_site_data'
                    parts = stripped.split()
                    if len(parts) >= len(cols):
                        rows.append(dict(zip(cols, parts[:len(cols)])))

            elif state == 'in_atom_site_data':
                if stripped == 'loop_' or stripped.startswith('data_'):
                    break
                parts = stripped.split()
                if len(parts) >= len(cols):
                    rows.append(dict(zip(cols, parts[:len(cols)])))

            elif state == 'in_other_loop':
                if stripped == 'loop_':
                    state = 'in_loop_start'
                    cols = []
                elif not stripped.startswith('_'):
                    # Data rows of a non-atom_site loop; keep consuming
                    pass

    return rows


# ============================================================================
# Binding site analysis
# ============================================================================

def compute_binding_sites(cif_path, distance_threshold=5.0):
    """
    Return (binding_site_Z_str, binding_site_Y_str) where each is a
    comma-separated string of residue entries "<seq_id>:<residue_name>"
    for protein residues with any atom within distance_threshold Å of
    the respective ligand.

    Returns ('', '') if the CIF file is missing or cannot be parsed.
    """
    if not os.path.exists(cif_path):
        return '', ''

    try:
        atoms = parse_atom_site(cif_path)
    except Exception:
        return '', ''

    if not atoms:
        return '', ''

    # Required columns — skip file if any are absent
    required = {'label_asym_id', 'label_seq_id', 'label_comp_id', 'Cartn_x', 'Cartn_y', 'Cartn_z'}
    if not required.issubset(atoms[0].keys()):
        return '', ''

    prot_atoms = []   # list of (seq_id, res_name, x, y, z)
    z_coords = []     # list of (x, y, z)
    y_coords = []     # list of (x, y, z)

    for a in atoms:
        chain = a['label_asym_id']
        try:
            x = float(a['Cartn_x'])
            y = float(a['Cartn_y'])
            z = float(a['Cartn_z'])
        except ValueError:
            continue

        if chain == 'A':
            seq_id = a['label_seq_id']
            res_name = a['label_comp_id']
            # Skip if residue number is not an integer (e.g. '.')
            try:
                seq_id = int(seq_id)
            except ValueError:
                continue
            prot_atoms.append((seq_id, res_name, x, y, z))
        elif chain == 'Z':
            z_coords.append((x, y, z))
        elif chain == 'Y':
            y_coords.append((x, y, z))

    if not prot_atoms:
        return '', ''

    prot_xyz = np.array([(a[2], a[3], a[4]) for a in prot_atoms], dtype=float)

    def residues_within(ligand_xyz_list):
        if not ligand_xyz_list:
            return ''
        lig_xyz = np.array(ligand_xyz_list, dtype=float)
        # prot_xyz: (N, 3), lig_xyz: (M, 3)
        diff = prot_xyz[:, np.newaxis, :] - lig_xyz[np.newaxis, :, :]  # (N, M, 3)
        dists = np.sqrt((diff ** 2).sum(axis=2))                        # (N, M)
        min_dist = dists.min(axis=1)                                    # (N,)
        close_residues = set()
        for idx, d in enumerate(min_dist):
            if d <= distance_threshold:
                seq_id, res_name = prot_atoms[idx][0], prot_atoms[idx][1]
                close_residues.add((seq_id, res_name))
        return ','.join(f"{sid}:{rname}" for sid, rname in sorted(close_residues))

    site_z = residues_within(z_coords)
    site_y = residues_within(y_coords)
    return site_z, site_y


# ============================================================================
# Score extraction
# ============================================================================

def _safe_val(v):
    """Return -1 if value is None, else the value itself."""
    return -1 if v is None else v


def read_confidences(summary_file, plddt_file):
    """
    Parse summary_confidences and confidences JSON files for a 3-chain run
    (A=protein, Z=ligand1, Y=ligand2).

    Returns a dict with per-ligand and global metrics, or None on failure.
    """
    try:
        with open(summary_file) as f:
            sc = json.load(f)
        with open(plddt_file) as f:
            conf = json.load(f)
    except Exception:
        return None

    chain_iptm = sc.get('chain_iptm', [])
    chain_ptm = sc.get('chain_ptm', [])
    pair_iptm = sc.get('chain_pair_iptm', [])
    pair_pae = sc.get('chain_pair_pae_min', [])

    atom_chain_ids = conf.get('atom_chain_ids', [])
    atom_plddts = conf.get('atom_plddts', [])

    def avg_plddt(chain_id):
        vals = [atom_plddts[i] for i, c in enumerate(atom_chain_ids) if c == chain_id]
        return float(np.mean(vals)) if vals else -1.0

    def ligand_metrics(idx):
        """Extract metrics for the ligand at chain index idx (1=Z, 2=Y)."""
        c_iptm = _safe_val(chain_iptm[idx]) if idx < len(chain_iptm) else -1
        c_ptm = _safe_val(chain_ptm[idx]) if idx < len(chain_ptm) else -1

        # pair_iptm[idx] is the row for this ligand; column 0 is vs protein A
        row_iptm = pair_iptm[idx] if idx < len(pair_iptm) else []
        p_iptm = _safe_val(row_iptm[0]) if row_iptm and row_iptm[0] is not None else -1

        row_pae = pair_pae[idx] if idx < len(pair_pae) else []
        p_pae = _safe_val(row_pae[0]) if row_pae and row_pae[0] is not None else -1

        return c_iptm, p_iptm, p_pae, c_ptm

    z_iptm, z_pair_iptm, z_pair_pae, z_ptm = ligand_metrics(1)
    y_iptm, y_pair_iptm, y_pair_pae, y_ptm = ligand_metrics(2)

    plddt_z = avg_plddt('Z')
    plddt_y = avg_plddt('Y')

    return {
        # Ligand Z
        'chain_iptm_Z': z_iptm,
        'pair_iptm_Z': z_pair_iptm,
        'pair_pae_min_Z': z_pair_pae,
        'chain_ptm_Z': z_ptm,
        'plddt_Z': plddt_z,
        # Ligand Y
        'chain_iptm_Y': y_iptm,
        'pair_iptm_Y': y_pair_iptm,
        'pair_pae_min_Y': y_pair_pae,
        'chain_ptm_Y': y_ptm,
        'plddt_Y': plddt_y,
        # Global
        'iptm': sc.get('iptm', -1),
        'ptm': sc.get('ptm', -1),
        'ranking_score': sc.get('ranking_score', -1),
    }


# ============================================================================
# Main scoring loop
# ============================================================================

def gen_scores_df(results_path, distance_threshold):
    columns = [
        'SMILES_Z', 'SMILES_Y', 'folder',
        'chain_iptm_Z', 'pair_iptm_Z', 'pair_pae_min_Z', 'chain_ptm_Z', 'plddt_Z',
        'chain_iptm_Y', 'pair_iptm_Y', 'pair_pae_min_Y', 'chain_ptm_Y', 'plddt_Y',
        'iptm', 'ptm', 'ranking_score',
        'binding_site_Z', 'binding_site_Y',
    ]

    records = []

    for folder in tqdm(sorted(os.listdir(results_path))):
        folder_path = os.path.join(results_path, folder)
        if not os.path.isdir(folder_path):
            continue
        if '_' in folder or 'token' in folder:
            continue

        data_file = os.path.join(folder_path, folder + '_data.json')
        summary_file = os.path.join(folder_path, folder + '_summary_confidences.json')
        plddt_file = os.path.join(folder_path, folder + '_confidences.json')
        cif_file = os.path.join(folder_path, folder + '_model.cif')

        if not os.path.exists(data_file):
            print(f"skip {folder}: missing _data.json")
            continue
        if not os.path.exists(summary_file) or not os.path.exists(plddt_file):
            print(f"skip {folder}: missing confidence files")
            continue

        # Read SMILES from the input data JSON
        try:
            with open(data_file) as f:
                data = json.load(f)
            smiles_z, smiles_y = '', ''
            for seq in data.get('sequences', []):
                if 'ligand' in seq:
                    lid = seq['ligand'].get('id', '')
                    if lid == 'Z':
                        smiles_z = seq['ligand'].get('smiles', '')
                    elif lid == 'Y':
                        smiles_y = seq['ligand'].get('smiles', '')
        except Exception:
            print(f"skip {folder}: failed to parse _data.json")
            continue

        metrics = read_confidences(summary_file, plddt_file)
        if metrics is None:
            print(f"skip {folder}: failed to parse confidence files")
            continue

        site_z, site_y = compute_binding_sites(cif_file, distance_threshold)

        row = [
            smiles_z, smiles_y, folder,
            metrics['chain_iptm_Z'], metrics['pair_iptm_Z'], metrics['pair_pae_min_Z'],
            metrics['chain_ptm_Z'], metrics['plddt_Z'],
            metrics['chain_iptm_Y'], metrics['pair_iptm_Y'], metrics['pair_pae_min_Y'],
            metrics['chain_ptm_Y'], metrics['plddt_Y'],
            metrics['iptm'], metrics['ptm'], metrics['ranking_score'],
            site_z, site_y,
        ]
        records.append(row)

    return pd.DataFrame(records, columns=columns)


def analyze(results_folder, output_dir, distance_threshold):
    df = gen_scores_df(results_folder, distance_threshold)
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, 'summary.csv')
    df.to_csv(output_path, index=False)
    print(f"Wrote {len(df)} rows to {output_path}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Summarize AF3 competitive folding scores (1 protein + 2 ligands)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('--af3-results-folder', required=True,
                        help="Folder containing AF3 competitive prediction outputs")
    parser.add_argument('--output-dir', required=True,
                        help="Directory to write summary.csv")
    parser.add_argument('--distance-threshold', type=float, default=5.0,
                        help="Ångström cutoff for binding site residue detection")
    args = parser.parse_args()

    analyze(args.af3_results_folder, args.output_dir, args.distance_threshold)
