# Import modules
import sys, platform
import json
from prody import *
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
import rdkit
from vina import Vina
import subprocess
import os
import pandas as pd
import argparse
from pdbfixer import PDBFixer
import openmm
from openmm.app import PDBFile, Modeller, ForceField
from openmm.app.simulation import Simulation
from openmm import LangevinIntegrator
from openmm.unit import kelvin, picosecond
from openmm import Platform

# Helper functions
def add_chain_identifier(pdb_file, chain_id='A'):
    """
    Check if PDB file has chain identifiers and add them if missing.
    Creates a new file with '_chain.pdb' suffix if modifications are needed.

    Args:
        pdb_file: Path to the PDB file
        chain_id: Chain identifier to add (default 'A')

    Returns:
        Path to the PDB file (original if already has chains, new file if modified)
    """
    with open(pdb_file, 'r') as f:
        lines = f.readlines()

    # Check if chain identifiers are present
    has_chain = False
    needs_chain = False

    for line in lines:
        if line.startswith(('ATOM', 'HETATM')):
            # Chain ID is at position 21 (0-indexed) in PDB format
            if len(line) > 21 and line[21].strip():
                has_chain = True
                break
            else:
                needs_chain = True

    if has_chain or not needs_chain:
        print(f"Chain identifiers already present in {pdb_file}")
        return pdb_file

    # Add chain identifiers
    print(f"Adding chain identifier '{chain_id}' to {pdb_file}")
    output_lines = []

    for line in lines:
        if line.startswith(('ATOM', 'HETATM')):
            # PDB format: positions 0-21 are fixed, position 21 is chain ID
            if len(line) > 21:
                # Insert chain ID at position 21
                modified_line = line[:21] + chain_id + line[22:]
                output_lines.append(modified_line)
            else:
                # Line is too short, pad and add chain ID
                padded_line = line.rstrip('\n').ljust(21) + chain_id + '\n'
                output_lines.append(padded_line)
        else:
            output_lines.append(line)

    # Create output filename
    base_name = pdb_file.rsplit('.pdb', 1)[0]
    output_file = f"{base_name}_chain.pdb"

    with open(output_file, 'w') as f:
        f.writelines(output_lines)

    print(f"Created {output_file} with chain identifiers")
    return output_file

def locate_file(from_path = None, query_path = None, query_name = "query file"):

    if not from_path or not query_path:
        raise ValueError("Must specify from_path and query_path")

    possible_path = list(from_path.glob(query_path))

    if not possible_path:
        raise FileNotFoundError(f"Cannot find {query_name} from {from_path} by {query_path}")

    return_which = (
        f"using {query_name} at:\n"
        f"{possible_path[0]}\n"
    )
    print(return_which)

    return possible_path[0]

def get_ca_center_of_mass(pdb_file):
    """
    Calculate center of mass of CA atoms from a PDB file.

    Args:
        pdb_file: Path to PDB file

    Returns:
        numpy array of shape (3,) with x, y, z coordinates in Angstroms
    """
    import numpy as np

    ca_coords = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                ca_coords.append([x, y, z])

    if not ca_coords:
        raise ValueError(f"No CA atoms found in {pdb_file}")

    ca_coords = np.array(ca_coords)
    com = np.mean(ca_coords, axis=0)
    return com

def get_ca_center_of_mass_from_openmm(topology, positions):
    """
    Calculate center of mass of CA atoms from OpenMM topology and positions.

    Args:
        topology: OpenMM Topology object
        positions: OpenMM positions (Quantity with units)

    Returns:
        numpy array of shape (3,) with x, y, z coordinates in Angstroms
    """
    import numpy as np
    from openmm.unit import nanometer

    # Convert positions to numpy array in Angstroms
    pos_array = np.array(positions.value_in_unit(nanometer)) * 10.0  # nm to Angstrom

    ca_indices = []
    for atom in topology.atoms():
        if atom.name == 'CA':
            ca_indices.append(atom.index)

    if not ca_indices:
        raise ValueError("No CA atoms found in topology")

    ca_coords = pos_array[ca_indices]
    com = np.mean(ca_coords, axis=0)
    return com

def energy_minimize_openmm(input_pdb, output_pdb, pH=7.0):
    """
    Energy minimization using OpenMM with coordinate preservation.

    Args:
        input_pdb: Path to input PDB file
        output_pdb: Path to output minimized PDB file
        pH: pH for adding hydrogens (default 7.0)

    Returns:
        Path to the minimized PDB file
    """
    import numpy as np
    from openmm.unit import nanometer, angstrom

    print(f"Energy minimizing {input_pdb} with OpenMM...")

    # Calculate original CA center of mass before PDBFixer
    original_ca_com = get_ca_center_of_mass(input_pdb)
    print(f"Original CA COM: {original_ca_com}")

    # 1) Load & fix
    fixer = PDBFixer(filename=input_pdb)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH)

    # 2) Build OpenMM system
    topology = fixer.topology
    positions = fixer.positions

    forcefield = ForceField("amber14-all.xml", "amber14/tip3p.xml")
    modeller = Modeller(topology, positions)
    # For receptor-only minimization, we typically skip adding solvent here.

    system = forcefield.createSystem(
        modeller.topology,
        #nonbondedMethod=openmm.NoCutoff
    )

    integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picosecond)

    platform = Platform.getPlatformByName("CPU")  # or "CUDA" / "OpenCL"
    simulation = Simulation(modeller.topology, system, integrator, platform)
    simulation.context.setPositions(modeller.positions)

    # 3) Minimize
    print("Minimizing...")
    simulation.minimizeEnergy()

    # 4) Get minimized positions
    state = simulation.context.getState(getPositions=True)
    min_positions = state.getPositions()

    # Calculate minimized CA center of mass
    minimized_ca_com = get_ca_center_of_mass_from_openmm(modeller.topology, min_positions)
    print(f"Minimized CA COM: {minimized_ca_com}")

    # Calculate translation vector to restore original position
    translation_vector = original_ca_com - minimized_ca_com
    translation_magnitude = np.linalg.norm(translation_vector)
    print(f"Translation vector: {translation_vector}")
    print(f"Translation magnitude: {translation_magnitude:.2f} Angstroms")

    # Warn if translation is very large
    if translation_magnitude > 100.0:
        print(f"WARNING: Large translation detected ({translation_magnitude:.2f} A). "
              f"Please verify results.")

    # Apply translation correction
    # Convert positions to numpy array in Angstroms, apply translation, convert back
    pos_array = np.array(min_positions.value_in_unit(nanometer)) * 10.0  # nm to Angstrom
    pos_array += translation_vector  # Apply translation
    corrected_positions = pos_array * angstrom  # Convert back to OpenMM Quantity

    # Write corrected PDB
    with open(output_pdb, "w") as f:
        PDBFile.writeFile(modeller.topology, corrected_positions, f)

    print(f"Written coordinate-corrected minimized structure to {output_pdb}")
    print(f"Applied translation correction: {translation_vector}")

    return output_pdb

def parse_boxes(args):
    """
    Parse box specifications from CLI args.
    Returns list of (center, size) tuples.
    """
    if args.boxes:
        vals = json.loads(args.boxes)
        # Normalize flat [6] to [[6]]
        if vals and not isinstance(vals[0], (list, tuple)):
            vals = [vals]
        boxes = []
        for v in vals:
            if len(v) != 6:
                raise ValueError(f"Each box must have 6 values [cx,cy,cz,sx,sy,sz], got {len(v)}")
            boxes.append((v[:3], v[3:]))
        return boxes
    elif args.box_center and args.box_size:
        return [(args.box_center, args.box_size)]
    else:
        raise ValueError("Must provide --boxes or both --box-center and --box-size")

def prepare_pdb(pdb_file, mk_prepare_receptor, outdir, centers, docking_box_size, box_idx=0):
    import numpy as np

    # Check and add chain identifiers if missing
    pdb_file = add_chain_identifier(pdb_file)

    reduce_opts = "approach=add\nadd_flip_movers=True"
    env = os.environ.copy()
    env["MMTBX_CCP4_MONOMER_LIB"] = geostd_path

    # Default name of reduce output...
    tmp_prefix = pdb_file.split(".pdb")[0].split("/")[-1]
    prepare_inPDB = f"{tmp_prefix}FH.pdb"
    output_filename = os.path.join(outdir, prepare_inPDB)

    # Compute original CA COM before reduce2
    original_ca_com = get_ca_center_of_mass(pdb_file)
    print(f"Original PDB CA COM: {original_ca_com}")

    # Run reduce2 once (cached by checking if output exists)
    if not os.path.exists(output_filename):
        print("Attempting meeko & reduce2 on original PDB...")
        subprocess.run(['python', reduce2, pdb_file, f"output.filename={output_filename}", reduce_opts],
                              env=env, capture_output=True)
    else:
        print(f"Reusing cached reduce2 output: {output_filename}")

    # Detect coordinate shift introduced by reduce2
    reduce2_ca_com = get_ca_center_of_mass(output_filename)
    coord_shift = reduce2_ca_com - original_ca_com
    shift_magnitude = np.linalg.norm(coord_shift)
    print(f"reduce2 CA COM: {reduce2_ca_com}")
    print(f"Coordinate shift from reduce2: {coord_shift} (magnitude: {shift_magnitude:.2f} A)")

    # Adjust box centers to follow the protein
    adjusted_centers = [
        centers[0] + coord_shift[0],
        centers[1] + coord_shift[1],
        centers[2] + coord_shift[2],
    ]
    if shift_magnitude > 1.0:
        print(f"Adjusting box centers: {list(centers)} -> {adjusted_centers}")

    # Size in each dimension
    center_x, center_y, center_z = adjusted_centers
    size_x, size_y, size_z = docking_box_size

    # mk_prepare_receptor runs per box (flexible residues are box-dependent)
    prepare_output = f"{outdir}/{tmp_prefix}FH_box{box_idx}"
    result = subprocess.run(["python", mk_prepare_receptor, "-i", output_filename, "-o", prepare_output, "-p", "-v", "--box_center", str(center_x), str(center_y), str(center_z), "--box_size", str(size_x), str(size_y), str(size_z)])

    # If reduce2 fails, try with energy minimization first
    if result.returncode != 0:
        print("Meeko failed on original PDB. Attempting energy minimization first...")
        minimized_pdb = os.path.join(outdir, f"{tmp_prefix}_minimized.pdb")
        if not os.path.exists(minimized_pdb):
            energy_minimize_openmm(pdb_file, minimized_pdb, pH=7.0)

        print("Retrying meeko & reduce2 on minimized PDB...")
        output_filename2 = os.path.join(outdir, f"{tmp_prefix}_minimizedFH.pdb")
        if not os.path.exists(output_filename2):
            subprocess.run(['python', reduce2, minimized_pdb, f"output.filename={output_filename2}", reduce_opts],
                                  env=env, capture_output=True)

        # Recompute shift for the fallback path
        minimized_ca_com = get_ca_center_of_mass(minimized_pdb)
        reduce2_min_ca_com = get_ca_center_of_mass(output_filename2)
        fallback_shift = reduce2_min_ca_com - minimized_ca_com
        # Total shift: original -> minimized is already corrected by energy_minimize_openmm,
        # so we only need the reduce2 shift on the minimized PDB
        # But minimized PDB was already corrected to match original coords,
        # so total shift from original box coords is just fallback_shift
        adjusted_centers = [
            centers[0] + fallback_shift[0],
            centers[1] + fallback_shift[1],
            centers[2] + fallback_shift[2],
        ]
        fallback_mag = np.linalg.norm(fallback_shift)
        print(f"Fallback reduce2 shift: {fallback_shift} (magnitude: {fallback_mag:.2f} A)")
        if fallback_mag > 1.0:
            print(f"Adjusting box centers (fallback): {list(centers)} -> {adjusted_centers}")

        center_x, center_y, center_z = adjusted_centers
        result = subprocess.run(["python", mk_prepare_receptor, "-i", output_filename2, "-o", prepare_output, "-p", "-v", "--box_center", str(center_x), str(center_y), str(center_z), "--box_size", str(size_x), str(size_y), str(size_z)])

        if result.returncode != 0:
            print("Meeko failed even after energy minimization!")
            print("STDERR:", result.stderr.decode())
            raise RuntimeError("Meeko failed to process the PDB file")
    else:
        print(f"Meeko & reduce2 succeeded for box{box_idx}.")

    # Return prefix and adjusted centers so caller uses corrected coordinates for docking
    return tmp_prefix, adjusted_centers

def prepare_ligand(lig_input, pH, scrub, mk_prepare_ligand, outdir):
    ligand_Smiles = lig_input

    args = ""
    skip_tautomer = True
    if skip_tautomer:
        args += "--skip_tautomer"
    skip_acidbase = False 
    if skip_acidbase:
        args += "--skip_acidbase"

    # Write scrubbed protomer(s) and conformer(s) to SDF
    ligandSDF = os.path.join(outdir, "prepared_ligand.sdf")
    pdbqt_out = os.path.join(outdir, "prepared_ligand.pdbqt")

    cmd = ["python", scrub, ligand_Smiles, "-o", ligandSDF,"--ph", str(pH),] + [args,]
    result = subprocess.run(cmd)
    
    subprocess.run(["python", mk_prepare_ligand, "-i", ligandSDF, "-o", pdbqt_out])
    return None

def is_already_docked(folder_path):
    """
    Check if docking has already been completed for this molecule.

    Args:
        folder_path: Path to the ligand folder

    Returns:
        True if docking output files exist, False otherwise
    """
    output_pdbqt = os.path.join(folder_path, 'docking_config_out.pdbqt')
    output_affinities = os.path.join(folder_path, 'docking_affinities.txt')

    # Check if both output files exist and are non-empty
    if os.path.exists(output_pdbqt) and os.path.exists(output_affinities):
        if os.path.getsize(output_pdbqt) > 0 and os.path.getsize(output_affinities) > 0:
            return True
    return False

def dock(receptorPDBQT, ligandPDBQT, centers, docking_box_size):

    v = Vina(sf_name='vina')
    v.set_receptor(rigid_pdbqt_filename=receptorPDBQT)
    v.set_ligand_from_file(ligandPDBQT)
    v.compute_vina_maps(center=centers, box_size=docking_box_size)

    # Dock the ligand
    v.dock(exhaustiveness=32, n_poses=5)
    v.write_poses('docking_config_out.pdbqt', n_poses=5, overwrite=True)
    affinities = v.energies()
    with open('docking_affinities.txt', 'w') as f:
        for i, affinity in enumerate(affinities):
            f.write(f'{affinity[0]:.5f} \n')
    return None

if __name__ == "__main__":
    # Command-line setup
    exe_path = "/home/ubuntu/miniconda3/envs/vina/bin/"
    reduce2_path = "/home/ubuntu/miniconda3/pkgs/cctbx-base-2024.2-py38hbbc5a03_0/lib/python3.8/site-packages/mmtbx/command_line/"
    geostd_p = "/home/ubuntu/Applications/"

    # Args
    parser = argparse.ArgumentParser(
        description="script to run autodock Vina.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--pdb", type=str, help="the protein pdb file (used when receptor pdbqt is not provided)")
    parser.add_argument("--receptor-pdbqt", type=str, help="prepared receptor pdbqt file")
    parser.add_argument("--smiles", type=str, help="the SMILES file")
    parser.add_argument("--boxes", type=str, help="JSON string of box specs: [[cx,cy,cz,sx,sy,sz], ...]")
    parser.add_argument("--box-center", type=float, nargs="+", help="(deprecated) center of the docking box")
    parser.add_argument("--box-size", type=float, nargs="+", help="(deprecated) the size of the docking box")
    parser.add_argument("--output", type=str, help="output folder")
    parser.add_argument("--smiles-col", type=str, default="SMILES", help="the name of the SMILES col")
    parser.add_argument("--skip-docked", action="store_true", help="skip molecules that are already docked")

    args = parser.parse_args()

    if not args.receptor_pdbqt and not args.pdb:
        parser.error("Must provide either --receptor-pdbqt or --pdb")

    def validate_root_block(pdbqt_file):
        root_line = None
        end_line = None

        with open(pdbqt_file, "r") as f:
            for lineno, line in enumerate(f, start=1):
                token = line.strip()
                if token == "ROOT" and root_line is None:
                    root_line = lineno
                elif token == "ENDROOT":
                    end_line = lineno

        if root_line is None or end_line is None:
            raise ValueError(f"Missing ROOT/ENDROOT markers in receptor PDBQT: {pdbqt_file}")
        if end_line <= root_line:
            raise ValueError(f"Invalid ROOT/ENDROOT ordering in receptor PDBQT: {pdbqt_file}")
        if end_line - root_line <= 1:
            raise ValueError(f"No lines between ROOT and ENDROOT in receptor PDBQT: {pdbqt_file}")

    scrub = locate_file(from_path=Path(exe_path), query_path="scrub.py", query_name="scrub.py")
    mk_prepare_ligand = locate_file(from_path=Path(exe_path), query_path="mk_prepare_ligand.py", query_name="mk_prepare_ligand.py")

    pH = 7.4
    boxes = parse_boxes(args)
    outdir = args.output
    os.makedirs(outdir, exist_ok=True)

    # Prepare receptor data for each box once (before ligand loop)
    receptor_pdbqts = {}
    adjusted_centers_per_box = {}

    if args.receptor_pdbqt:
        receptor_pdbqt = os.path.abspath(args.receptor_pdbqt)
        if not os.path.exists(receptor_pdbqt):
            raise FileNotFoundError(f"Receptor PDBQT not found: {receptor_pdbqt}")

        validate_root_block(receptor_pdbqt)
        print(f"Using pre-prepared receptor PDBQT: {receptor_pdbqt}")

        for box_idx, (center, _size) in enumerate(boxes):
            receptor_pdbqts[box_idx] = receptor_pdbqt
            adjusted_centers_per_box[box_idx] = center
    else:
        pdb_file = args.pdb
        mk_prepare_receptor = locate_file(
            from_path=Path(exe_path),
            query_path="mk_prepare_receptor.py",
            query_name="mk_prepare_receptor.py",
        )
        reduce2 = locate_file(from_path=Path(reduce2_path), query_path="reduce2.py", query_name="reduce2.py")
        geostd_path = locate_file(from_path=Path(geostd_p), query_path="geostd", query_name="geostd")

        for box_idx, (center, size) in enumerate(boxes):
            tmp_prefix, adj_center = prepare_pdb(
                pdb_file,
                mk_prepare_receptor,
                outdir,
                center,
                size,
                box_idx=box_idx,
            )
            receptor_pdbqts[box_idx] = os.path.join(outdir, f"{tmp_prefix}FH_box{box_idx}.pdbqt")
            adjusted_centers_per_box[box_idx] = adj_center

    parent_folder = args.output
    tmp_df = pd.read_csv(args.smiles)

    # Track success and failures per (ligand, box) pair
    successful_dockings = []
    failed_dockings = []
    skipped_dockings = []

    for i, lig in enumerate(tmp_df[args.smiles_col]):
        lig_folder = os.path.join(parent_folder, f"lig{i}")
        os.makedirs(lig_folder, exist_ok=True)

        # Prepare ligand once per molecule
        ligand_prepared = False
        pdbqt_out = os.path.join(lig_folder, "prepared_ligand.pdbqt")

        for box_idx, (center, size) in enumerate(boxes):
            box_folder = os.path.join(lig_folder, f"box{box_idx}")
            os.makedirs(box_folder, exist_ok=True)

            # Check if already docked (only if skip-docked flag is set)
            if args.skip_docked and is_already_docked(box_folder):
                print(f"Skipping lig{i}/box{box_idx}: Already docked")
                skipped_dockings.append((i, box_idx))
                continue

            try:
                # Prepare ligand once (cached across boxes)
                if not ligand_prepared:
                    print(f"Preparing lig{i}...")
                    prepare_ligand(lig, pH, scrub, mk_prepare_ligand, lig_folder)
                    ligand_prepared = True

                os.chdir(box_folder)
                receptorPDBQT_path = receptor_pdbqts[box_idx]

                print(f"Docking lig{i} to box{box_idx}...")
                dock(receptorPDBQT_path, pdbqt_out, adjusted_centers_per_box[box_idx], size)

                successful_dockings.append((i, box_idx))
                print(f"Successfully completed docking for lig{i}/box{box_idx}")

            except Exception as e:
                failed_dockings.append((i, box_idx))
                print(f"ERROR: Failed to dock lig{i}/box{box_idx}: {str(e)}")
                print("Continuing with next...")
                error_log = os.path.join(box_folder, "docking_error.log")
                with open(error_log, "w") as f:
                    f.write(f"Error during docking: {str(e)}\n")

    # Print summary
    total_pairs = len(tmp_df) * len(boxes)
    print("\n" + "=" * 60)
    print("DOCKING SUMMARY")
    print("=" * 60)
    print(f"Total molecules: {len(tmp_df)}")
    print(f"Boxes per molecule: {len(boxes)}")
    print(f"Total (ligand, box) pairs: {total_pairs}")
    print(f"Successfully docked: {len(successful_dockings)}")
    if args.skip_docked:
        print(f"Skipped (already docked): {len(skipped_dockings)}")
    print(f"Failed: {len(failed_dockings)}")

    if failed_dockings:
        print(f"\nFailed pairs: {failed_dockings}")
    if args.skip_docked and skipped_dockings:
        print(f"Skipped pairs: {skipped_dockings}")
    print("=" * 60)
