#!/usr/bin/env python3
"""Prepare, validate, and cache shared ColabFold A3M files."""

import argparse
import csv
import gzip
import hashlib
import json
import shutil
from pathlib import Path


MANIFEST_FORMAT_VERSION = 2


def chain_ids(count):
    letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    result = []
    for index in range(count):
        if index < len(letters):
            result.append(letters[index])
        else:
            result.append(letters[index // len(letters) - 1] + letters[index % len(letters)])
    return result


def load_chains(sequences_csv, protein, name_col="protein_name", seq_col="sequence"):
    with Path(sequences_csv).open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        if name_col not in fieldnames:
            for candidate in ("protein_name", "name"):
                if candidate in fieldnames:
                    name_col = candidate
                    break
        if name_col not in fieldnames or seq_col not in fieldnames:
            raise ValueError(
                f"{sequences_csv} must contain a protein-name column and {seq_col!r}"
            )
        row = next((item for item in reader if str(item[name_col]) == protein), None)
    if row is None:
        raise ValueError(f"Protein {protein!r} not found in {sequences_csv}")
    if row.get(seq_col) is None:
        raise ValueError(
            f"Protein {protein!r} has no {seq_col!r} value in {sequences_csv}"
        )
    sequences = [part.strip().upper() for part in str(row[seq_col]).split(",")]
    sequences = [sequence for sequence in sequences if sequence]
    if not sequences:
        raise ValueError(f"Protein {protein!r} has no non-empty chains")
    return list(zip(chain_ids(len(sequences)), sequences))


def sha256_text(value):
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def parse_a3m(path):
    records = []
    description = None
    sequence_parts = []
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            # ColabFold may prefix an A3M with multimer length/cardinality
            # metadata even when each query was searched independently.
            if line.startswith("#") and description is None:
                continue
            if line.startswith(">"):
                if description is not None:
                    records.append((description, "".join(sequence_parts)))
                description = line[1:]
                sequence_parts = []
            else:
                if description is None:
                    raise ValueError(f"{path} has sequence data before its first header")
                sequence_parts.append(line)
    if description is not None:
        records.append((description, "".join(sequence_parts)))
    if not records:
        raise ValueError(f"{path} contains no A3M records")
    return records


def validate_a3m(path, query):
    records = parse_a3m(path)
    first = "".join(character for character in records[0][1] if not character.islower())
    if first.upper() != query:
        raise ValueError(f"{path} first sequence does not match the query")
    for description, sequence in records:
        aligned = "".join(character for character in sequence if not character.islower())
        if len(aligned) != len(query):
            raise ValueError(
                f"{path} record {description!r} has {len(aligned)} aligned columns; "
                f"expected {len(query)}"
            )
    return len(records)


def sequence_digest(chains):
    canonical = "\n".join(f"{chain_id}:{sequence}" for chain_id, sequence in chains)
    return sha256_text(canonical)


def command_prepare(args):
    chains = load_chains(args.sequences_csv, args.protein, args.name_col, args.seq_col)
    work_dir = Path(args.work_dir).resolve()
    input_dir = work_dir / "input"
    search_dir = work_dir / "search"
    input_dir.mkdir(parents=True, exist_ok=True)
    search_dir.mkdir(parents=True, exist_ok=True)

    representatives = {}
    search_rows = []
    chain_rows = []
    for chain_id, sequence in chains:
        digest = sha256_text(sequence)
        representative = representatives.get(digest)
        if representative is None:
            representative = chain_id
            representatives[digest] = representative
            fasta_path = input_dir / f"{args.protein}_{representative}.fasta"
            fasta_path.write_text(
                f">{args.protein}_{representative}\n{sequence}\n", encoding="utf-8"
            )
            output_dir = search_dir / representative
            search_rows.append(
                {
                    "representative_chain": representative,
                    "sequence_sha256": digest,
                    "fasta_path": str(fasta_path),
                    "output_dir": str(output_dir),
                }
            )
        chain_rows.append(
            {
                "chain_id": chain_id,
                "sequence": sequence,
                "sequence_sha256": digest,
                "representative_chain": representative,
            }
        )

    plan = {
        "protein": args.protein,
        "sequences_csv": str(Path(args.sequences_csv).resolve()),
        "sequence_digest": sequence_digest(chains),
        "chains": chain_rows,
        "searches": search_rows,
    }
    plan_path = work_dir / "search_plan.json"
    plan_path.write_text(json.dumps(plan, indent=2) + "\n", encoding="utf-8")
    with (work_dir / "search_plan.tsv").open("w", encoding="utf-8") as handle:
        for row in search_rows:
            handle.write(
                f"{row['representative_chain']}\t{row['fasta_path']}\t{row['output_dir']}\n"
            )
    print(plan_path)


def find_a3m(search_dir):
    candidates = sorted(
        list(Path(search_dir).rglob("*.a3m"))
        + list(Path(search_dir).rglob("*.a3m.gz"))
    )
    if not candidates:
        raise FileNotFoundError(f"No A3M output found under {search_dir}")
    preferred = [path for path in candidates if "coverage" not in path.name]
    return (preferred or candidates)[0]


def command_finalize(args):
    work_dir = Path(args.work_dir).resolve()
    msa_dir = Path(args.msa_dir).resolve()
    plan = json.loads((work_dir / "search_plan.json").read_text(encoding="utf-8"))
    normalized_dir = msa_dir / "colabfold"
    normalized_dir.mkdir(parents=True, exist_ok=True)

    source_by_representative = {}
    for search in plan["searches"]:
        source_by_representative[search["representative_chain"]] = find_a3m(
            search["output_dir"]
        )

    manifest_chains = []
    for chain in plan["chains"]:
        source = source_by_representative[chain["representative_chain"]]
        record_count = validate_a3m(source, chain["sequence"])
        # OpenFold3 uses an A3M file's parent directory as its representative
        # identifier, while Boltz requires identical protein sequences to use
        # the exact same MSA path. Use one directory/path per unique sequence:
        # duplicate chains share it, and different sequences remain distinct.
        representative_dir = (
            normalized_dir
            / f"{plan['protein']}_{chain['representative_chain']}"
        )
        representative_dir.mkdir(parents=True, exist_ok=True)
        destination = representative_dir / "colabfold_main.a3m"
        temporary = destination.with_suffix(destination.suffix + ".tmp")
        if source.suffix == ".gz":
            with gzip.open(source, "rt", encoding="utf-8") as source_handle:
                temporary.write_text(source_handle.read(), encoding="utf-8")
        else:
            shutil.copyfile(source, temporary)
        validate_a3m(temporary, chain["sequence"])
        temporary.replace(destination)
        manifest_chains.append(
            {
                "chain_id": chain["chain_id"],
                "sequence_sha256": chain["sequence_sha256"],
                "representative_chain": chain["representative_chain"],
                "msa_path": str(destination),
                "record_count": record_count,
            }
        )

    manifest = {
        "format_version": MANIFEST_FORMAT_VERSION,
        "source": "colabfold_public_server",
        "protein": plan["protein"],
        "sequences_csv": plan["sequences_csv"],
        "sequence_digest": plan["sequence_digest"],
        "colabfold_command": args.colabfold_command,
        "colabfold_version": args.colabfold_version,
        "chains": manifest_chains,
    }
    manifest_path = msa_dir / "manifest.json"
    temporary_manifest = manifest_path.with_suffix(".json.tmp")
    temporary_manifest.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    temporary_manifest.replace(manifest_path)
    (msa_dir / "msa.done").write_text(f"{manifest['sequence_digest']}\n", encoding="utf-8")
    print(manifest_path)


def check_cache(msa_dir, sequences_csv, protein, name_col, seq_col):
    msa_dir = Path(msa_dir).resolve()
    manifest_path = msa_dir / "manifest.json"
    token_path = msa_dir / "msa.done"
    if not manifest_path.is_file() or not token_path.is_file():
        return False
    chains = load_chains(sequences_csv, protein, name_col, seq_col)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("format_version") != MANIFEST_FORMAT_VERSION:
        return False
    if manifest.get("protein") != protein:
        return False
    if manifest.get("sequence_digest") != sequence_digest(chains):
        return False
    if token_path.read_text(encoding="utf-8").strip() != manifest["sequence_digest"]:
        return False
    entries = {entry["chain_id"]: entry for entry in manifest.get("chains", [])}
    if set(entries) != {chain_id for chain_id, _ in chains}:
        return False
    for chain_id, sequence in chains:
        entry = entries[chain_id]
        if entry.get("sequence_sha256") != sha256_text(sequence):
            return False
        path = Path(entry.get("msa_path", ""))
        if not path.is_file():
            return False
        validate_a3m(path, sequence)
    return True


def command_check(args):
    if not check_cache(
        args.msa_dir,
        args.sequences_csv,
        args.protein,
        args.name_col,
        args.seq_col,
    ):
        raise SystemExit(1)
    print(Path(args.msa_dir).resolve() / "manifest.json")


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    for name in ("prepare", "check"):
        subparser = subparsers.add_parser(name)
        subparser.add_argument("--sequences-csv", required=True)
        subparser.add_argument("--protein", required=True)
        subparser.add_argument("--name-col", default="protein_name")
        subparser.add_argument("--seq-col", default="sequence")
        if name == "prepare":
            subparser.add_argument("--work-dir", required=True)
            subparser.set_defaults(func=command_prepare)
        else:
            subparser.add_argument("--msa-dir", required=True)
            subparser.set_defaults(func=command_check)
    finalize = subparsers.add_parser("finalize")
    finalize.add_argument("--work-dir", required=True)
    finalize.add_argument("--msa-dir", required=True)
    finalize.add_argument("--colabfold-command", required=True)
    finalize.add_argument("--colabfold-version", default="unknown")
    finalize.set_defaults(func=command_finalize)
    return parser


def main():
    args = build_parser().parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
