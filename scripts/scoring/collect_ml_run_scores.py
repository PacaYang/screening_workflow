#!/usr/bin/env python3
"""Backfill Nesso, Protenix, and OpenFold3 scores across a benchmark run."""

import argparse
import csv
from pathlib import Path

from ml_score_common import (
    COMMON_COLUMNS, METHODS, SCORE_COLUMNS, atomic_write_csv, collect_method,
)


def read_rows(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", required=True, type=Path)
    args = parser.parse_args()
    run_root = args.run_root
    predictions_path = run_root / "predictions.csv"
    if not predictions_path.is_file():
        parser.error(f"predictions.csv not found: {predictions_path}")
    predictions = read_rows(predictions_path)
    if not predictions:
        parser.error("predictions.csv is empty")
    context_rows = {}
    for prediction in predictions:
        context = prediction["context"]
        context_rows.setdefault(context, []).append(prediction)

    summaries = {}
    for context in sorted(context_rows):
        context_dir = run_root / context
        selected = context_dir / "initial_screening" / "selected.csv"
        if not selected.is_file():
            raise ValueError(f"selected.csv not found: {selected}")
        for method in METHODS:
            method_dir = context_dir / "fine_screening" / method
            summary_path, counts = collect_method(method, method_dir / "output", selected, method_dir)
            rows = read_rows(summary_path)
            if len(rows) != len(context_rows[context]):
                raise ValueError(f"row count mismatch in {summary_path}")
            summaries[(context, method)] = rows
            print(f"{context}/{method}: {counts}")

    columns = list(predictions[0])
    method_columns = {}
    for method in METHODS:
        prefix = method.lower() + "_"
        known = set(COMMON_COLUMNS) | set(SCORE_COLUMNS[method])
        extra = {field for context in context_rows for row in summaries[(context, method)]
                 for field in row if field not in known}
        fields = list(COMMON_COLUMNS) + list(SCORE_COLUMNS[method]) + sorted(extra)
        method_columns[method] = [field for field in fields if field not in {"compound_index", "SMILES"}]
        columns.extend(prefix + field for field in method_columns[method])

    merged = []
    merged_by_id = {}
    for context, context_predictions in context_rows.items():
        per_method = {}
        for method in METHODS:
            rows = summaries[(context, method)]
            indexed = {int(row["compound_index"]): row for row in rows}
            if len(indexed) != len(rows):
                raise ValueError(f"duplicate compound index in {context}/{method}")
            per_method[method] = indexed
        for index, prediction in enumerate(context_predictions):
            expected_id = f"{context}__{index}"
            if prediction["prediction_id"] != expected_id:
                raise ValueError(f"prediction order does not match selected.csv: {expected_id}")
            row = dict(prediction)
            for method in METHODS:
                scores = per_method[method][index]
                if scores["SMILES"] != prediction["SMILES"]:
                    raise ValueError(f"SMILES mismatch for {expected_id}/{method}")
                for field in method_columns[method]:
                    row[method.lower() + "_" + field] = scores.get(field, "")
            merged_by_id[expected_id] = row
    for prediction in predictions:
        merged.append(merged_by_id[prediction["prediction_id"]])
    output = run_root / "compiled_scores.csv"
    atomic_write_csv(output, merged, columns)
    print(f"Wrote {len(merged)} predictions to {output}")


if __name__ == "__main__":
    main()
