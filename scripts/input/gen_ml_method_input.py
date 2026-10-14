#!/usr/bin/env python3
"""Generate canonical protein/ligand inputs for AF3-like structure models."""

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from seq_utils import chain_ids, ligand_ids, parse_chains


def load_target(sequences_csv, protein, name_col, seq_col):
    frame = pd.read_csv(sequences_csv)
    if name_col not in frame.columns:
        for candidate in ("protein_name", "name"):
            if candidate in frame.columns:
                name_col = candidate
                break
    if name_col not in frame.columns or seq_col not in frame.columns:
        raise ValueError(f"{sequences_csv} must contain {name_col!r} and {seq_col!r}")
    rows = frame.loc[frame[name_col].astype(str) == protein]
    if rows.empty:
        raise ValueError(f"Protein {protein!r} not found in {sequences_csv}")
    chains = parse_chains(rows.iloc[0][seq_col])
    if not chains:
        raise ValueError(f"Protein {protein!r} has no non-empty sequences")
    return chains


def load_msa_paths(manifest_path, protein, chains):
    if not manifest_path:
        return None
    manifest_file = Path(manifest_path).resolve()
    manifest = json.loads(manifest_file.read_text(encoding="utf-8"))
    if manifest.get("protein") != protein:
        raise ValueError(f"MSA manifest is for {manifest.get('protein')!r}, not {protein!r}")
    entries = {entry["chain_id"]: entry for entry in manifest.get("chains", [])}
    paths = []
    for chain_id, sequence in zip(chain_ids(len(chains)), chains):
        entry = entries.get(chain_id)
        if not entry:
            raise ValueError(f"MSA manifest has no entry for chain {chain_id}")
        digest = hashlib.sha256(sequence.upper().encode("utf-8")).hexdigest()
        if entry.get("sequence_sha256") != digest:
            raise ValueError(f"MSA manifest sequence hash is stale for chain {chain_id}")
        msa_path = Path(entry.get("msa_path", ""))
        if not msa_path.is_file():
            raise ValueError(f"MSA file not found for chain {chain_id}: {msa_path}")
        paths.append(str(msa_path.resolve()))
    return paths


def make_input(
    method, protein, chains, ligands, mode, index, template_dir=None, msa_paths=None
):
    proteins = []
    for chain_id, sequence in zip(chain_ids(len(chains)), chains):
        entry = {
            "id": chain_id,
            "sequence": sequence,
            "pairedMsa": "",
            "templates": [],
        }
        if msa_paths:
            entry["unpairedMsaPath"] = msa_paths[len(proteins)]
        else:
            entry["unpairedMsa"] = ""
        if template_dir:
            cif = Path(template_dir) / f"{protein}_{chain_id}.cif"
            if cif.exists():
                entry["templates"] = [{
                    "mmcifPath": str(cif.resolve()),
                    "queryIndices": list(range(len(sequence))),
                    "templateIndices": list(range(len(sequence))),
                }]
        proteins.append({"protein": entry})

    payload = {
        "name": f"{protein}_{index}" if mode != "prefold" else protein,
        "dialect": "alphafold3",
        "version": 2,
        "modelSeeds": [99],
        "sequences": proteins,
        "screening": {
            "method": method,
            "protein": protein,
            "compound_index": index,
            "mode": mode,
        },
    }
    # AlphaFold 3 rejects unknown top-level keys. Other adapters may use the
    # screening metadata internally, but AF3 inputs must stay on its schema.
    if method.lower() == "af3":
        payload.pop("screening")
    if ligands:
        payload["sequences"].extend(
            {"ligand": {"id": lid, "smiles": smiles}}
            for lid, smiles in zip(ligand_ids(len(ligands)), ligands)
        )
    return payload


def make_nesso_input(protein, chains, ligands, mode, index):
    """Build NESSO's version-1 YAML schema from the canonical AF3 payload."""
    sequences = [
        {"protein": {"id": chain_id, "sequence": sequence}}
        for chain_id, sequence in zip(chain_ids(len(chains)), chains)
    ]
    if ligands:
        sequences.extend(
            {"ligand": {"id": ligand_id, "smiles": smiles}}
            for ligand_id, smiles in zip(ligand_ids(len(ligands)), ligands)
        )
    payload = {
        "version": 1,
        "name": f"{protein}_{index}" if mode != "prefold" else protein,
        "sequences": sequences,
    }
    if ligands:
        payload["properties"] = [{"affinity": {"binder": ligand_ids(1)[0]}}]
    return payload


def make_protenix_input(protein, chains, ligands, mode, index, msa_paths=None):
    """Build Protenix's entity schema from the canonical input."""
    sequences = []
    for chain_index, sequence in enumerate(chains):
        protein_chain = {"sequence": sequence, "count": 1}
        if msa_paths:
            protein_chain["unpairedMsaPath"] = msa_paths[chain_index]
            protein_chain["pairedMsa"] = ""
        sequences.append({"proteinChain": protein_chain})
    if ligands:
        sequences.extend(
            {"ligand": {"ligand": smiles, "count": 1}}
            for smiles in ligands
        )
    return [{
        "name": f"{protein}_{index}" if mode != "prefold" else protein,
        "sequences": sequences,
    }]


def make_openfold3_input(protein, chains, ligands, mode, index, msa_paths=None):
    """Build OpenFold3's named query-set schema from the canonical input."""
    query = {
        "chains": [
            {
                "molecule_type": "protein",
                "chain_ids": [chain_id],
                "sequence": sequence,
                **(
                    {"main_msa_file_paths": [msa_paths[chain_index]]}
                    if msa_paths
                    else {}
                ),
            }
            for chain_index, (chain_id, sequence) in enumerate(
                zip(chain_ids(len(chains)), chains)
            )
        ],
        "use_msas": bool(msa_paths),
        "use_paired_msas": False,
        "use_main_msas": bool(msa_paths),
    }
    query["chains"].extend(
        {
            "molecule_type": "ligand",
            "chain_ids": [ligand_id],
            "smiles": smiles,
        }
        for ligand_id, smiles in zip(ligand_ids(len(ligands)), ligands)
    )
    return {
        "seeds": [99],
        "queries": {
            f"{protein}_{index}" if mode != "prefold" else protein: query
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--method", required=True)
    parser.add_argument("--mode", choices=("prefold", "complex", "template"), required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--input-csv", required=True)
    parser.add_argument("--protein-name", required=True)
    parser.add_argument("--smiles-file")
    parser.add_argument("--smiles-col", default="SMILES")
    parser.add_argument("--name-col", default="protein_name")
    parser.add_argument("--seq-col", default="sequence")
    parser.add_argument("--template-dir")
    parser.add_argument("--msa-manifest")
    parser.add_argument("--n-ligands", type=int, default=1)
    args = parser.parse_args()

    chains = load_target(args.input_csv, args.protein_name, args.name_col, args.seq_col)
    msa_paths = load_msa_paths(args.msa_manifest, args.protein_name, chains)
    if args.mode == "prefold":
        compound_rows = [(0, [])]
    else:
        if not args.smiles_file:
            raise ValueError("--smiles-file is required outside prefold mode")
        compounds = pd.read_csv(args.smiles_file)
        if args.smiles_col not in compounds.columns:
            raise ValueError(f"{args.smiles_file} has no {args.smiles_col!r} column")
        smiles = ["" if pd.isna(value) else str(value).strip() for value in compounds[args.smiles_col].tolist()]
        if any(not value for value in smiles):
            raise ValueError("SMILES input contains an empty value")
        group_size = max(1, args.n_ligands)
        compound_rows = [
            (start, smiles[start:start + group_size])
            for start in range(0, len(smiles), group_size)
        ]

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "input_manifest.jsonl").open("w", encoding="utf-8") as manifest:
        for chunk_index, (index, ligands) in enumerate(compound_rows):
            payload = make_input(
                args.method, args.protein_name, chains, ligands, args.mode, index,
                args.template_dir if args.mode == "template" else None,
                msa_paths,
            )
            filename = f"{args.protein_name}.json" if args.mode == "prefold" else f"{args.protein_name}_{chunk_index}.json"
            path = output_dir / filename
            path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
            input_path = path
            if args.method.lower() == "nesso":
                input_path = path.with_suffix(".yaml")
                input_path.write_text(
                    json.dumps(
                        make_nesso_input(
                            args.protein_name, chains, ligands, args.mode, index
                        ),
                        indent=2,
                    )
                    + "\n",
                    encoding="utf-8",
                )
            elif args.method.lower() == "protenix":
                input_path = path.with_suffix(".protenix.json")
                input_path.write_text(
                    json.dumps(
                        make_protenix_input(
                            args.protein_name, chains, ligands, args.mode, index,
                            msa_paths,
                        ),
                        indent=2,
                    )
                    + "\n",
                    encoding="utf-8",
                )
            elif args.method.lower() == "openfold3":
                input_path = path.with_suffix(".openfold.json")
                input_path.write_text(
                    json.dumps(
                        make_openfold3_input(
                            args.protein_name, chains, ligands, args.mode, index,
                            msa_paths,
                        ),
                        indent=2,
                    )
                    + "\n",
                    encoding="utf-8",
                )
            manifest.write(json.dumps({
                "method": args.method,
                "protein": args.protein_name,
                "compound_index": index,
                "smiles": ligands[0] if ligands else "",
                "ligands": ligands,
                "mode": args.mode,
                "input": str(input_path),
            }) + "\n")
    print(f"Wrote {len(compound_rows)} {args.method} {args.mode} input(s) to {output_dir}")


if __name__ == "__main__":
    main()
