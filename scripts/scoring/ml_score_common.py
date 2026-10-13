#!/usr/bin/env python3
"""Collect native Nesso, Protenix, and OpenFold3 scores by selected-row index."""

from __future__ import annotations

import csv
import json
import math
import os
import re
import tempfile
from pathlib import Path


METHODS = ("Nesso", "Protenix", "OpenFold3")
COMMON_COLUMNS = (
    "compound_index", "folder", "SMILES", "status", "reason",
    "seed", "sample", "sample_count", "score_json",
)
SCORE_COLUMNS = {
    "Nesso": (
        "affinity_pred_value", "affinity_pred_value1", "affinity_pred_value2",
        "affinity_logits_binary", "affinity_probability_binary",
        "entropy_pp", "entropy_pl", "entropy_ll",
        "entropy_crop_pp", "entropy_crop_pl", "entropy_crop_ll",
    ),
    "Protenix": (
        "plddt", "gpde", "ptm", "iptm", "has_clash", "disorder",
        "ranking_score", "num_recycles", "ligand_chain_ptm",
        "ligand_chain_iptm", "ligand_chain_plddt", "ligand_chain_gpde",
        "ligand_protein_pair_iptm_max", "ligand_protein_pair_iptm_global_max",
        "ligand_protein_pair_plddt_max", "ligand_protein_pair_gpde_min",
    ),
    "OpenFold3": (
        "avg_plddt", "gpde", "iptm", "ptm", "disorder", "has_clash",
        "sample_ranking_score", "ligand_chain_ptm",
        "ligand_protein_pair_iptm_max",
        "ligand_protein_bespoke_iptm_max",
    ),
}
RANK_KEY = {
    "Nesso": "affinity_pred_value",
    "Protenix": "ranking_score",
    "OpenFold3": "sample_ranking_score",
}
SAMPLE_RE = re.compile(r"_sample_(\d+)(?:_confidences_aggregated)?\.json$")
BATCH_RE = re.compile(r"^batch_(\d+)$")
JOB_RE = re.compile(r"^job_(\d+)$")
PAIR_RE = re.compile(r"^\(([^,]+),\s*([^)]+)\)$")


def read_selected(path: Path) -> list[str]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or "SMILES" not in reader.fieldnames:
            raise ValueError(f"selected CSV has no SMILES column: {path}")
        smiles = [(row.get("SMILES") or "").strip() for row in reader]
    if any(not value for value in smiles):
        raise ValueError(f"selected CSV contains an empty SMILES value: {path}")
    return smiles


def atomic_write_csv(path: Path, rows: list[dict], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            "w", newline="", encoding="utf-8", dir=path.parent,
            prefix=f".{path.name}.", suffix=".tmp", delete=False,
        ) as handle:
            temporary = Path(handle.name)
            writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, path)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def read_input_manifest(path: Path, expected_count: int) -> dict[int, tuple[int, int, str]]:
    """Map selected-row index to (job index, group size, declared SMILES)."""
    if not path.exists():
        return {i: (i, 1, "") for i in range(expected_count)}
    mapping = {}
    with path.open(encoding="utf-8") as handle:
        for job_index, line in enumerate(handle):
            if not line.strip():
                continue
            record = json.loads(line)
            start = int(record["compound_index"])
            ligands = record.get("ligands", [])
            if not isinstance(ligands, list) or not ligands:
                raise ValueError(f"invalid ligands in {path}, record {job_index}")
            for offset, ligand in enumerate(ligands):
                index = start + offset
                if index >= expected_count or index in mapping:
                    raise ValueError(f"invalid compound index {index} in {path}")
                mapping[index] = (job_index, len(ligands), str(ligand))
    if set(mapping) != set(range(expected_count)):
        raise ValueError(f"input manifest does not cover selected CSV: {path}")
    return mapping


def read_job_results(output_dir: Path, job_count: int) -> dict[int, str]:
    """Return latest terminal batch status for each job index."""
    ledger = output_dir / "job_results.tsv"
    job_ids = output_dir / "job_ids.txt"
    if not ledger.exists() or not job_ids.exists() or job_count == 0:
        return {}
    if not job_ids.read_text(encoding="utf-8").strip():
        return {}
    batch_status = {}
    batch_sizes = []
    with ledger.open(encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            match = BATCH_RE.fullmatch(parts[3])
            if match and parts[4] in {"success", "failed"}:
                batch_status[int(match.group(1))] = parts[4]
                if len(parts) >= 7 and parts[6].isdigit():
                    batch_sizes.append(int(parts[6]))
    if not batch_sizes:
        return {}
    batch_size = max(batch_sizes)
    if batch_size <= 0:
        return {}
    return {index: batch_status[index // batch_size]
            for index in range(job_count) if index // batch_size in batch_status}


def finite_number(value: object) -> float | None:
    if isinstance(value, bool):
        return None
    if not isinstance(value, (int, float)):
        return None
    number = float(value)
    return number if math.isfinite(number) else None


def scalar_scores(data: dict) -> dict:
    return {key: value for key, value in data.items()
            if isinstance(value, bool) or finite_number(value) is not None}


def indexed_value(values: object, index: int) -> object:
    if isinstance(values, list) and index < len(values):
        value = values[index]
        if finite_number(value) is not None:
            return value
    return ""


def matrix_extreme(matrix: object, index: int, mode: str) -> object:
    if not isinstance(matrix, list) or index >= len(matrix):
        return ""
    row = matrix[index]
    if not isinstance(row, list):
        return ""
    values = [value for i, value in enumerate(row)
              if i != index and finite_number(value) is not None]
    if not values:
        return ""
    return max(values) if mode == "max" else min(values)


def openfold_pair_extreme(pairs: object, ligand: str) -> object:
    if not isinstance(pairs, dict):
        return ""
    values = []
    for key, value in pairs.items():
        match = PAIR_RE.fullmatch(key)
        if match and ligand in match.groups() and finite_number(value) is not None:
            values.append(value)
    return max(values) if values else ""


def ligand_scores(method: str, data: dict) -> dict:
    if method == "Protenix":
        chains = data.get("chain_ptm")
        if not isinstance(chains, list) or len(chains) < 2:
            return {}
        index = len(chains) - 1  # Generated inputs put the sole ligand last.
        return {
            "ligand_chain_ptm": indexed_value(chains, index),
            "ligand_chain_iptm": indexed_value(data.get("chain_iptm"), index),
            "ligand_chain_plddt": indexed_value(data.get("chain_plddt"), index),
            "ligand_chain_gpde": indexed_value(data.get("chain_gpde"), index),
            "ligand_protein_pair_iptm_max": matrix_extreme(data.get("chain_pair_iptm"), index, "max"),
            "ligand_protein_pair_iptm_global_max": matrix_extreme(data.get("chain_pair_iptm_global"), index, "max"),
            "ligand_protein_pair_plddt_max": matrix_extreme(data.get("chain_pair_plddt"), index, "max"),
            "ligand_protein_pair_gpde_min": matrix_extreme(data.get("chain_pair_gpde"), index, "min"),
        }
    if method == "OpenFold3":
        chains = data.get("chain_ptm")
        if not isinstance(chains, dict) or len(chains) < 2:
            return {}
        ligand = "Z" if "Z" in chains else list(chains)[-1]
        return {
            "ligand_chain_ptm": chains.get(ligand, ""),
            "ligand_protein_pair_iptm_max": openfold_pair_extreme(data.get("chain_pair_iptm"), ligand),
            "ligand_protein_bespoke_iptm_max": openfold_pair_extreme(data.get("bespoke_iptm"), ligand),
        }
    return {}


def candidate_files(method: str, job_dir: Path, query: str) -> list[Path]:
    if method == "Nesso":
        pattern = "predictions/*/affinity.json"
    elif method == "Protenix":
        pattern = "*/seed_*/predictions/*_summary_confidence_sample_*.json"
    else:
        pattern = "*/seed_*/*_confidences_aggregated.json"
    files = sorted(job_dir.glob(pattern))
    if method == "Nesso":
        return [path for path in files if path.parent.name == query]
    if method == "Protenix":
        return [path for path in files if path.name.startswith(query + "_summary_confidence_sample_")]
    return [path for path in files if path.name.startswith(query + "_seed_")]


def sample_metadata(method: str, path: Path) -> tuple[int | str, int | str]:
    if method == "Nesso":
        return "", ""
    seed_match = re.search(r"seed_(\d+)", str(path.parent))
    sample_match = SAMPLE_RE.search(path.name)
    if not seed_match or not sample_match:
        raise ValueError(f"unrecognized seed/sample in {path}")
    return int(seed_match.group(1)), int(sample_match.group(1))


def collect_method(method: str, results_folder: Path, selected_csv: Path,
                   output_dir: Path) -> tuple[Path, dict[str, int]]:
    if method not in METHODS:
        raise ValueError(f"unsupported method: {method}")
    smiles = read_selected(selected_csv)
    context = selected_csv.parent.parent.name
    manifest = results_folder.parent / "input" / "input_manifest.jsonl"
    mapping = read_input_manifest(manifest, len(smiles))
    job_count = max((entry[0] for entry in mapping.values()), default=-1) + 1
    terminal = read_job_results(results_folder, job_count)
    rows = []
    counts = {}
    extra_scalar_keys = set()
    for index, selected_smiles in enumerate(smiles):
        job_index, group_size, declared_smiles = mapping[index]
        folder = f"job_{job_index}"
        job_dir = results_folder / folder
        row = {"compound_index": index, "folder": folder, "SMILES": selected_smiles,
               "status": "pending", "reason": "", "seed": "", "sample": "",
               "sample_count": 0, "score_json": ""}
        if group_size != 1:
            row.update(status="unsupported_group", reason="multiple ligands in one model input")
        elif declared_smiles and declared_smiles != selected_smiles:
            row.update(status="input_mismatch", reason="manifest SMILES differs from selected.csv")
        elif job_dir.is_dir() and (job_dir / "smiles.txt").exists() and \
                (job_dir / "smiles.txt").read_text(encoding="utf-8").strip() != selected_smiles:
            row.update(status="input_mismatch", reason="job smiles.txt differs from selected.csv")
        else:
            query = f"{context}_{job_index}"
            paths = candidate_files(method, job_dir, query) if job_dir.is_dir() else []
            valid = []
            problems = []
            for path in paths:
                try:
                    data = json.loads(path.read_text(encoding="utf-8"))
                    if not isinstance(data, dict):
                        raise ValueError("JSON root is not an object")
                    primary = finite_number(data.get(RANK_KEY[method]))
                    if primary is None:
                        raise ValueError(f"missing or invalid {RANK_KEY[method]}")
                    if method == "Nesso" and finite_number(data.get("affinity_probability_binary")) is None:
                        raise ValueError("missing or invalid affinity_probability_binary")
                    seed, sample = sample_metadata(method, path)
                    valid.append((primary, seed, sample, path, data))
                except (OSError, ValueError, json.JSONDecodeError) as error:
                    problems.append(f"{path.name}: {error}")
            if valid:
                if method == "Nesso":
                    chosen = valid[0]
                else:
                    chosen = min(valid, key=lambda item: (-item[0], item[1], item[2], str(item[3])))
                _, seed, sample, path, data = chosen
                scores = scalar_scores(data)
                extra_scalar_keys.update(scores)
                row.update(scores)
                row.update(ligand_scores(method, data))
                row.update(status="scored", reason=(f"ignored {len(problems)} invalid sample(s)" if problems else ""),
                           seed=seed, sample=sample, sample_count=len(valid),
                           score_json=str(path.relative_to(results_folder)))
            elif job_index in terminal:
                if problems:
                    row.update(status="invalid_score", reason="; ".join(problems[:2]))
                elif terminal[job_index] == "failed":
                    row.update(status="failed_job", reason="job ledger reports failed")
                else:
                    row.update(status="missing_score", reason="terminal job has no native score JSON")
            else:
                row.update(status="pending", reason="awaiting terminal job result")
        rows.append(row)
        counts[row["status"]] = counts.get(row["status"], 0) + 1
    columns = list(COMMON_COLUMNS) + list(SCORE_COLUMNS[method])
    columns += sorted(extra_scalar_keys - set(columns))
    output_path = output_dir / "summary.csv"
    atomic_write_csv(output_path, rows, columns)
    return output_path, counts
