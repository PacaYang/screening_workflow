###################################
# Convert the sdf file from diffdock to mol file
# Used for PBSA MD step
###################################

from rdkit import Chem
from rdkit.Chem import AllChem
import os
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert SDF file to mol file")

    parser.add_argument(
        "-input",
        type=str,
        required=True,
        help="SDF file path"
    )
    args = parser.parse_args()

    sdf_file = args.input

    suppl = Chem.SDMolSupplier(sdf_file)
    for mol in suppl:
        if mol is None:
            continue

        frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)

        # Choose the largest fragment (most atoms)
        if not frags:
            continue

        largest_frag = max(frags, key=lambda m: m.GetNumAtoms())

        # Sanitize the fragment
        Chem.SanitizeMol(largest_frag)

        Chem.MolToMolFile(largest_frag, "ligand.mol")
