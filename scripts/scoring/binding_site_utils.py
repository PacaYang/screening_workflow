"""
Shared utilities for binding-site extraction from structure prediction outputs.

Supports:
  - mmCIF files (Boltz2, AF3): parse_atom_site + compute_binding_site_cif
  - PDB files (RoseTTAFold):   parse_pdb_atoms + compute_binding_site_pdb

Binding site definition: protein residues with any atom within `cutoff_ang`
Angstroms of any ligand atom (default 10.0 Å = 1 nm).

Output format: comma-separated "<resnum>:<resname>" sorted by residue number,
e.g. "123:LEU,124:ALA,130:GLY".  Empty string on failure or no residues found.
"""

import os
import numpy as np


# ============================================================================
# mmCIF parser
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
                else:
                    state = 'searching'

            elif state == 'in_atom_site_header':
                if stripped.startswith('_atom_site.'):
                    cols.append(stripped[len('_atom_site.'):])
                elif stripped == 'loop_':
                    break
                elif stripped.startswith('data_'):
                    break
                else:
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
                    pass

    return rows


def compute_binding_site_cif(cif_path, ligand_chain, protein_chain='A', cutoff_ang=10.0):
    """
    Return a comma-separated string of protein residues within `cutoff_ang`
    Angstroms of the ligand in `ligand_chain`.

    Format: "<resnum>:<resname>,..." sorted numerically by residue number.
    Returns '' if the file is missing, unparseable, or no residues are found.

    Parameters
    ----------
    cif_path : str
        Path to the mmCIF structure file.
    ligand_chain : str
        Chain ID of the ligand (e.g. 'Z' for Boltz2/AF3).
    protein_chain : str
        Chain ID of the protein (default 'A').
    cutoff_ang : float
        Distance cutoff in Angstroms (default 10.0 = 1 nm).
    """
    if not os.path.exists(cif_path):
        return ''

    try:
        atoms = parse_atom_site(cif_path)
    except Exception:
        return ''

    if not atoms:
        return ''

    required = {'label_asym_id', 'label_seq_id', 'label_comp_id', 'Cartn_x', 'Cartn_y', 'Cartn_z'}
    if not required.issubset(atoms[0].keys()):
        return ''

    prot_atoms = []   # (seq_id_int, res_name, x, y, z)
    lig_coords = []   # (x, y, z)

    for a in atoms:
        chain = a['label_asym_id']
        try:
            x = float(a['Cartn_x'])
            y = float(a['Cartn_y'])
            z = float(a['Cartn_z'])
        except ValueError:
            continue

        if chain == protein_chain:
            try:
                seq_id = int(a['label_seq_id'])
            except ValueError:
                continue
            prot_atoms.append((seq_id, a['label_comp_id'], x, y, z))
        elif chain == ligand_chain:
            lig_coords.append((x, y, z))

    if not prot_atoms or not lig_coords:
        return ''

    prot_xyz = np.array([(a[2], a[3], a[4]) for a in prot_atoms], dtype=float)
    lig_xyz = np.array(lig_coords, dtype=float)

    diff = prot_xyz[:, np.newaxis, :] - lig_xyz[np.newaxis, :, :]  # (N, M, 3)
    min_dist = np.sqrt((diff ** 2).sum(axis=2)).min(axis=1)         # (N,)

    close = set()
    for idx, d in enumerate(min_dist):
        if d <= cutoff_ang:
            seq_id, res_name = prot_atoms[idx][0], prot_atoms[idx][1]
            close.add((seq_id, res_name))

    return ','.join(f"{sid}:{rname}" for sid, rname in sorted(close))


# ============================================================================
# PDB parser
# ============================================================================

def parse_pdb_atoms(pdb_path):
    """
    Parse ATOM and HETATM records from a PDB file.

    Returns a list of dicts with keys:
        chain_id, res_seq (int), res_name, x, y, z
    Lines that cannot be parsed are skipped silently.
    """
    atoms = []
    with open(pdb_path) as f:
        for line in f:
            rec = line[:6].strip()
            if rec not in ('ATOM', 'HETATM'):
                continue
            try:
                chain_id = line[21]
                res_seq = int(line[22:26].strip())
                res_name = line[17:20].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except (ValueError, IndexError):
                continue
            atoms.append({
                'chain_id': chain_id,
                'res_seq': res_seq,
                'res_name': res_name,
                'x': x,
                'y': y,
                'z': z,
            })
    return atoms


def compute_binding_site_pdb(pdb_path, ligand_chain, protein_chain='A', cutoff_ang=10.0):
    """
    Return a comma-separated string of protein residues within `cutoff_ang`
    Angstroms of the ligand in `ligand_chain`, parsed from a PDB file.

    Format: "<resnum>:<resname>,..." sorted numerically by residue number.
    Returns '' if the file is missing, unparseable, or no residues are found.

    Parameters
    ----------
    pdb_path : str
        Path to the PDB structure file.
    ligand_chain : str
        Chain ID of the ligand (e.g. 'B' for RoseTTAFold).
    protein_chain : str
        Chain ID of the protein (default 'A').
    cutoff_ang : float
        Distance cutoff in Angstroms (default 10.0 = 1 nm).
    """
    if not os.path.exists(pdb_path):
        return ''

    try:
        atoms = parse_pdb_atoms(pdb_path)
    except Exception:
        return ''

    if not atoms:
        return ''

    prot_atoms = []   # (res_seq_int, res_name, x, y, z)
    lig_coords = []   # (x, y, z)

    for a in atoms:
        if a['chain_id'] == protein_chain:
            prot_atoms.append((a['res_seq'], a['res_name'], a['x'], a['y'], a['z']))
        elif a['chain_id'] == ligand_chain:
            lig_coords.append((a['x'], a['y'], a['z']))

    if not prot_atoms or not lig_coords:
        return ''

    prot_xyz = np.array([(a[2], a[3], a[4]) for a in prot_atoms], dtype=float)
    lig_xyz = np.array(lig_coords, dtype=float)

    diff = prot_xyz[:, np.newaxis, :] - lig_xyz[np.newaxis, :, :]  # (N, M, 3)
    min_dist = np.sqrt((diff ** 2).sum(axis=2)).min(axis=1)         # (N,)

    close = set()
    for idx, d in enumerate(min_dist):
        if d <= cutoff_ang:
            res_seq, res_name = prot_atoms[idx][0], prot_atoms[idx][1]
            close.add((res_seq, res_name))

    return ','.join(f"{sid}:{rname}" for sid, rname in sorted(close))
