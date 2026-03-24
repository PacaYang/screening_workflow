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
import shutil
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

def locate_obabel(exe_path):
    candidate = Path(exe_path) / "obabel"
    if candidate.exists():
        print(f"using obabel at:\n{candidate}\n")
        return str(candidate)

    which_path = shutil.which("obabel")
    if which_path:
        print(f"using obabel at:\n{which_path}\n")
        return which_path

    raise FileNotFoundError(
        "Cannot find obabel. Please install Open Babel in the vina_test environment."
    )

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
    return parse_boxes_optional(args, required=True)

def parse_boxes_optional(args, required=True):
    if args.boxes:
        vals = json.loads(args.boxes)
        # Normalize flat [6] to [[6]]
        if vals and not isinstance(vals[0], (list, tuple)):
            vals = [vals]
        boxes = []
        for v in vals:
            if len(v) != 6:
                raise ValueError(f"Each box must have 6 values [cx,cy,cz,sx,sy,sz], got {len(v)}")
            boxes.append((list(map(float, v[:3])), list(map(float, v[3:]))))
        return boxes
    if args.box_center and args.box_size:
        if len(args.box_center) != 3 or len(args.box_size) != 3:
            raise ValueError("--box-center and --box-size must each contain exactly 3 values")
        return [(list(map(float, args.box_center)), list(map(float, args.box_size)))]
    if required:
        raise ValueError("Must provide --boxes or both --box-center and --box-size")
    return None

def boxes_to_specs(boxes):
    specs = []
    for center, size in boxes:
        specs.append([
            float(center[0]), float(center[1]), float(center[2]),
            float(size[0]), float(size[1]), float(size[2]),
        ])
    return specs

def specs_to_boxes(specs):
    boxes = []
    for idx, spec in enumerate(specs):
        if len(spec) != 6:
            raise ValueError(f"Invalid box spec at index {idx}: expected 6 values, got {len(spec)}")
        boxes.append((
            [float(spec[0]), float(spec[1]), float(spec[2])],
            [float(spec[3]), float(spec[4]), float(spec[5])],
        ))
    return boxes

def extract_root_atom_rows(input_pdbqt, output_pdbqt):
    with open(input_pdbqt, "r") as f:
        lines = f.readlines()

    saw_root = False
    saw_endroot = False
    in_root = False
    root_atom_rows = []

    for line in lines:
        token = line.strip()
        if not saw_root and token == "ROOT":
            saw_root = True
            in_root = True
            continue
        if in_root and token == "ENDROOT":
            saw_endroot = True
            break
        if in_root and line.startswith(("ATOM", "HETATM")):
            root_atom_rows.append(line)

    if not saw_root:
        raise RuntimeError(f"ROOT marker not found in {input_pdbqt}")
    if not saw_endroot:
        raise RuntimeError(f"ENDROOT marker not found in {input_pdbqt}")
    if not root_atom_rows:
        raise RuntimeError(f"No ATOM/HETATM rows found between ROOT and ENDROOT in {input_pdbqt}")

    with open(output_pdbqt, "w") as f:
        f.writelines(root_atom_rows)

    print(
        f"Wrote ROOT-only receptor rows to {output_pdbqt} "
        f"({len(root_atom_rows)} ATOM/HETATM lines)"
    )

def prepare_receptor_pdbqt_with_obabel(pdb_file, obabel_exe, outdir):
    pdb_file = add_chain_identifier(pdb_file)
    tmp_prefix = os.path.basename(pdb_file).rsplit(".pdb", 1)[0]

    raw_pdbqt = os.path.join(outdir, f"{tmp_prefix}_obabel_raw.pdbqt")

    if os.path.exists(raw_pdbqt) and os.path.getsize(raw_pdbqt) > 0:
        print(f"Reusing cached receptor PDBQT: {raw_pdbqt}")
        return os.path.abspath(raw_pdbqt), tmp_prefix

    cmd = [obabel_exe, "-ipdb", pdb_file, "-opdbqt", "-O", raw_pdbqt, "-xr", "--addpolarh", "--partialcharge", "gasteiger"]
    print(f"Running Open Babel receptor conversion: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            "Open Babel receptor conversion failed.\n"
            f"STDOUT: {result.stdout}\nSTDERR: {result.stderr}"
        )

    if not os.path.exists(raw_pdbqt):
        raise RuntimeError(f"Open Babel did not produce expected PDBQT file: {raw_pdbqt}")

    atom_rows = 0
    with open(raw_pdbqt, "r") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                atom_rows += 1

    if atom_rows == 0:
        raise RuntimeError(f"Open Babel receptor PDBQT has no ATOM/HETATM rows: {raw_pdbqt}")

    print(f"Prepared full receptor PDBQT: {raw_pdbqt} ({atom_rows} ATOM/HETATM lines)")
    return os.path.abspath(raw_pdbqt), tmp_prefix

def prepare_receptors_for_boxes(pdb_file, obabel_exe, outdir, boxes):
    receptor_pdbqt_paths = {}
    adjusted_centers_per_box = {}
    receptor_pdbqt_path, _ = prepare_receptor_pdbqt_with_obabel(pdb_file, obabel_exe, outdir)

    for box_idx, (center, _) in enumerate(boxes):
        receptor_pdbqt_paths[box_idx] = receptor_pdbqt_path
        adjusted_centers_per_box[box_idx] = list(map(float, center))

    return receptor_pdbqt_paths, adjusted_centers_per_box

def write_receptor_manifest(prepared_dir, source_pdb, boxes, adjusted_centers_per_box, receptor_pdbqt_paths):
    manifest_path = os.path.join(prepared_dir, "receptor_manifest.json")

    boxes_input = boxes_to_specs(boxes)
    boxes_adjusted = []
    receptor_pdbqts = []

    for idx, (_, size) in enumerate(boxes):
        if idx not in adjusted_centers_per_box:
            raise ValueError(f"Missing adjusted center for box{idx}")
        if idx not in receptor_pdbqt_paths:
            raise ValueError(f"Missing receptor PDBQT path for box{idx}")

        adj = adjusted_centers_per_box[idx]
        boxes_adjusted.append([
            float(adj[0]), float(adj[1]), float(adj[2]),
            float(size[0]), float(size[1]), float(size[2]),
        ])
        receptor_pdbqts.append(os.path.basename(receptor_pdbqt_paths[idx]))

    manifest = {
        "source_pdb": os.path.abspath(source_pdb),
        "boxes_input": boxes_input,
        "boxes_adjusted": boxes_adjusted,
        "receptor_pdbqts": receptor_pdbqts,
        "receptor_prep_mode": "obabel_raw_full",
    }

    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)

    print(f"Wrote receptor manifest to {manifest_path}")
    return manifest_path

def _specs_match(a, b, tol=1e-3):
    if len(a) != len(b):
        return False
    for i in range(len(a)):
        if len(a[i]) != len(b[i]):
            return False
        for j in range(len(a[i])):
            if abs(float(a[i][j]) - float(b[i][j])) > tol:
                return False
    return True

def load_prepared_receptors(prepared_dir, requested_boxes=None):
    manifest_path = os.path.join(prepared_dir, "receptor_manifest.json")
    if not os.path.exists(manifest_path):
        raise FileNotFoundError(f"Prepared receptor manifest not found: {manifest_path}")

    with open(manifest_path, "r") as f:
        manifest = json.load(f)

    if "boxes_input" not in manifest or "boxes_adjusted" not in manifest or "receptor_pdbqts" not in manifest:
        raise ValueError(f"Invalid receptor manifest: missing required keys in {manifest_path}")

    boxes_input = specs_to_boxes(manifest["boxes_input"])
    boxes_adjusted = specs_to_boxes(manifest["boxes_adjusted"])

    receptor_entries = manifest["receptor_pdbqts"]
    if isinstance(receptor_entries, dict):
        ordered_keys = sorted(receptor_entries.keys(), key=lambda x: int(x))
        receptor_names = [receptor_entries[k] for k in ordered_keys]
    else:
        receptor_names = receptor_entries

    if not (len(boxes_input) == len(boxes_adjusted) == len(receptor_names)):
        raise ValueError(
            "Prepared receptor manifest has inconsistent box/receptor counts: "
            f"{len(boxes_input)} input, {len(boxes_adjusted)} adjusted, {len(receptor_names)} receptor files"
        )

    if requested_boxes is not None:
        requested_specs = boxes_to_specs(requested_boxes)
        manifest_specs = boxes_to_specs(boxes_input)
        if not _specs_match(requested_specs, manifest_specs):
            raise ValueError(
                "Requested boxes do not match prepared receptor manifest. "
                "Regenerate prepared receptors or pass matching --boxes."
            )
        boxes_to_use = requested_boxes
    else:
        boxes_to_use = boxes_input

    adjusted_centers_per_box = {}
    receptor_pdbqt_paths = {}

    for idx, ((adj_center, _), receptor_name) in enumerate(zip(boxes_adjusted, receptor_names)):
        receptor_path = receptor_name
        if not os.path.isabs(receptor_path):
            receptor_path = os.path.join(prepared_dir, receptor_path)
        if not os.path.exists(receptor_path):
            raise FileNotFoundError(f"Prepared receptor file missing for box{idx}: {receptor_path}")

        adjusted_centers_per_box[idx] = list(map(float, adj_center))
        receptor_pdbqt_paths[idx] = receptor_path

    print(f"Loaded prepared receptors from {prepared_dir}")
    return boxes_to_use, adjusted_centers_per_box, receptor_pdbqt_paths

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
    prepare_input_pdb = output_filename

    # Compute original CA COM before reduce2
    original_ca_com = get_ca_center_of_mass(pdb_file)
    print(f"Original PDB CA COM: {original_ca_com}")

    # Run reduce2 once (cached by checking if output exists)
    if not os.path.exists(output_filename):
        print("Attempting meeko & reduce2 on original PDB...")
        reduce_result = subprocess.run(['python', reduce2, pdb_file, f"output.filename={output_filename}", reduce_opts],
                              env=env, capture_output=True)
        if reduce_result.returncode != 0:
            print(f"reduce2 returned non-zero exit code: {reduce_result.returncode}")
            print("reduce2 STDERR:", reduce_result.stderr.decode() if reduce_result.stderr else "(no stderr)")
        if not os.path.exists(output_filename):
            print(f"WARNING: reduce2 did not produce expected output file: {output_filename}")
            prepare_input_pdb = pdb_file
    else:
        print(f"Reusing cached reduce2 output: {output_filename}")

    # Detect coordinate shift introduced by reduce2
    if os.path.exists(output_filename):
        reduce2_ca_com = get_ca_center_of_mass(output_filename)
        coord_shift = reduce2_ca_com - original_ca_com
        shift_magnitude = np.linalg.norm(coord_shift)
        print(f"reduce2 CA COM: {reduce2_ca_com}")
        print(f"Coordinate shift from reduce2: {coord_shift} (magnitude: {shift_magnitude:.2f} A)")
    else:
        coord_shift = np.array([0.0, 0.0, 0.0])
        shift_magnitude = 0.0
        print("Skipping reduce2 coordinate-shift correction because reduce2 output is unavailable.")

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
    mk_cmd = ["python", mk_prepare_receptor, "-i", prepare_input_pdb, "-o", prepare_output, "-p", "-v", "--box_center", str(center_x), str(center_y), str(center_z), "--box_size", str(size_x), str(size_y), str(size_z)]
    result = subprocess.run(mk_cmd, capture_output=True)

    # If mk_prepare_receptor fails, check for histidine ambiguity first
    if result.returncode != 0:
        import re
        stderr_text = result.stderr.decode() if result.stderr else ""
        print(f"mk_prepare_receptor failed (exit code {result.returncode}):")
        print("STDERR:", stderr_text)
        his_pattern = r"for residue_key='(.*?)', .* tied for fewest missing H: (\w+)"
        his_matches = re.findall(his_pattern, stderr_text)

        if his_matches:
            # Resolve ambiguous histidines by picking the first tied variant
            assignments = []
            for res_key, first_variant in his_matches:
                print(f"Resolving ambiguous histidine: {res_key} -> {first_variant}")
                assignments.append(f"{res_key}={first_variant}")
            set_template_arg = ",".join(assignments)
            print(f"Retrying mk_prepare_receptor with --set_template {set_template_arg}")
            result = subprocess.run(mk_cmd + ["--set_template", set_template_arg], capture_output=True)
            if result.returncode != 0:
                print("mk_prepare_receptor with --set_template failed:")
                print("STDERR:", result.stderr.decode() if result.stderr else "(no stderr)")

        if result.returncode != 0:
            # Try --allow_bad_res as a safety net before energy minimization
            print("Trying mk_prepare_receptor with --allow_bad_res...")
            allow_bad_cmd = mk_cmd + ["--allow_bad_res"]
            if his_matches:
                allow_bad_cmd += ["--set_template", set_template_arg]
            result = subprocess.run(allow_bad_cmd, capture_output=True)
            if result.returncode != 0:
                print("mk_prepare_receptor with --allow_bad_res failed:")
                print("STDERR:", result.stderr.decode() if result.stderr else "(no stderr)")

        if result.returncode != 0:
            print("Meeko failed on original PDB. Attempting energy minimization first...")
            minimized_pdb = os.path.join(outdir, f"{tmp_prefix}_minimized.pdb")
            if not os.path.exists(minimized_pdb):
                energy_minimize_openmm(pdb_file, minimized_pdb, pH=7.0)

            # Recompute shift for the fallback path
            minimized_ca_com = get_ca_center_of_mass(minimized_pdb)
            fallback_shift = minimized_ca_com - original_ca_com
            fallback_mag = np.linalg.norm(fallback_shift)
            print(f"Fallback shift (original -> minimized): {fallback_shift} (magnitude: {fallback_mag:.2f} A)")
            adjusted_centers = [
                centers[0] + fallback_shift[0],
                centers[1] + fallback_shift[1],
                centers[2] + fallback_shift[2],
            ]
            if fallback_mag > 1.0:
                print(f"Adjusting box centers (fallback): {list(centers)} -> {adjusted_centers}")
            center_x, center_y, center_z = adjusted_centers

            # Skip reduce2 on minimized PDB — OpenMM already adds hydrogens via addMissingHydrogens()
            print("Running mk_prepare_receptor directly on minimized PDB (already has hydrogens)...")
            mk_cmd2 = ["python", mk_prepare_receptor, "-i", minimized_pdb, "-o", prepare_output, "-p", "-v", "--box_center", str(center_x), str(center_y), str(center_z), "--box_size", str(size_x), str(size_y), str(size_z)]
            result = subprocess.run(mk_cmd2, capture_output=True)

            # Try histidine fix on minimized PDB too
            if result.returncode != 0:
                stderr_text2 = result.stderr.decode() if result.stderr else ""
                print("mk_prepare_receptor on minimized PDB failed:")
                print("STDERR:", stderr_text2)
                his_matches2 = re.findall(his_pattern, stderr_text2)
                if his_matches2:
                    assignments2 = []
                    for res_key, first_variant in his_matches2:
                        print(f"Resolving ambiguous histidine (minimized): {res_key} -> {first_variant}")
                        assignments2.append(f"{res_key}={first_variant}")
                    set_template_arg2 = ",".join(assignments2)
                    print(f"Retrying with --set_template {set_template_arg2}")
                    result = subprocess.run(mk_cmd2 + ["--set_template", set_template_arg2], capture_output=True)

            # Last resort: --allow_bad_res on minimized PDB
            if result.returncode != 0:
                print("Trying --allow_bad_res on minimized PDB...")
                allow_bad_cmd2 = mk_cmd2 + ["--allow_bad_res"]
                if his_matches2:
                    allow_bad_cmd2 += ["--set_template", set_template_arg2]
                result = subprocess.run(allow_bad_cmd2, capture_output=True)

            if result.returncode != 0:
                print("Meeko failed even after energy minimization!")
                print("STDERR:", result.stderr.decode() if result.stderr else "(no stderr)")
                print("STDOUT:", result.stdout.decode() if result.stdout else "(no stdout)")
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
    # Commandline scripts
    exe_path = "/home/yangl_pacagen_com/miniconda3/envs/vina_test/bin/"

    # Args
    parser = argparse.ArgumentParser(description="script to run autodock Vina.", \
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--pdb", type=str, help="the protein pdb file")
    parser.add_argument("--smiles", type=str, help="the SMILES file")
    parser.add_argument("--boxes", type=str, help="JSON string of box specs: [[cx,cy,cz,sx,sy,sz], ...]")
    parser.add_argument("--box-center", type=float, nargs="+",  help="(deprecated) center of the docking box")
    parser.add_argument("--box-size", type=float, nargs="+",  help="(deprecated) the size of the docking box")
    parser.add_argument('--output', type=str, help="output folder")
    parser.add_argument("--smiles-col", type=str, default='SMILES', help="the name of the SMILES col")
    parser.add_argument("--skip-docked", action='store_true', help="skip molecules that are already docked")
    parser.add_argument(
        "--prepare-receptor-only",
        action="store_true",
        help="prepare receptor assets/manifest and exit without docking ligands",
    )
    parser.add_argument(
        "--prepared-receptor-dir",
        type=str,
        help="directory containing receptor_manifest.json and prepared receptor PDBQT files",
    )

    args = parser.parse_args()

    if args.prepare_receptor_only and args.prepared_receptor_dir:
        parser.error("--prepare-receptor-only cannot be used with --prepared-receptor-dir")

    needs_obabel = args.prepare_receptor_only or not args.prepared_receptor_dir
    needs_ligand_tools = not args.prepare_receptor_only

    obabel_exe = None
    scrub = None
    mk_prepare_ligand = None

    if needs_obabel:
        obabel_exe = locate_obabel(exe_path)
    if needs_ligand_tools:
        scrub = locate_file(from_path = Path(exe_path), query_path = "scrub.py", query_name = "scrub.py")
        mk_prepare_ligand = locate_file(from_path = Path(exe_path), query_path = "mk_prepare_ligand.py", query_name = "mk_prepare_ligand.py")

    needs_boxes = args.prepare_receptor_only or not args.prepared_receptor_dir
    boxes = parse_boxes_optional(args, required=needs_boxes)

    if args.prepare_receptor_only:
        if not args.pdb:
            parser.error("--pdb is required with --prepare-receptor-only")
        if not args.output:
            parser.error("--output is required with --prepare-receptor-only")

        os.makedirs(args.output, exist_ok=True)
        receptor_pdbqt_paths, adjusted_centers_per_box = prepare_receptors_for_boxes(
            args.pdb, obabel_exe, args.output, boxes
        )
        write_receptor_manifest(
            args.output, args.pdb, boxes, adjusted_centers_per_box, receptor_pdbqt_paths
        )
        print("Receptor preprocessing completed.")
        sys.exit(0)

    if not args.smiles:
        parser.error("--smiles is required for docking mode")
    if not args.output:
        parser.error("--output is required for docking mode")

    parent_folder = args.output
    os.makedirs(parent_folder, exist_ok=True)
    pH = 7.4

    if args.prepared_receptor_dir:
        boxes, adjusted_centers_per_box, receptor_pdbqt_paths = load_prepared_receptors(
            args.prepared_receptor_dir, requested_boxes=boxes
        )
    else:
        if not args.pdb:
            parser.error("--pdb is required when --prepared-receptor-dir is not provided")
        receptor_pdbqt_paths, adjusted_centers_per_box = prepare_receptors_for_boxes(
            args.pdb, obabel_exe, parent_folder, boxes
        )
        write_receptor_manifest(
            parent_folder, args.pdb, boxes, adjusted_centers_per_box, receptor_pdbqt_paths
        )

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
        pdbqt_out = os.path.join(lig_folder, 'prepared_ligand.pdbqt')

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
                receptorPDBQT_path = receptor_pdbqt_paths[box_idx]

                print(f"Docking lig{i} to box{box_idx}...")
                dock(receptorPDBQT_path, pdbqt_out, adjusted_centers_per_box[box_idx], size)

                successful_dockings.append((i, box_idx))
                print(f"Successfully completed docking for lig{i}/box{box_idx}")

            except Exception as e:
                failed_dockings.append((i, box_idx))
                print(f"ERROR: Failed to dock lig{i}/box{box_idx}: {str(e)}")
                print(f"Continuing with next...")
                error_log = os.path.join(box_folder, 'docking_error.log')
                with open(error_log, 'w') as f:
                    f.write(f"Error during docking: {str(e)}\n")

    # Print summary
    total_pairs = len(tmp_df) * len(boxes)
    print("\n" + "="*60)
    print("DOCKING SUMMARY")
    print("="*60)
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
    print("="*60)
