#!/usr/bin/env python3
"""Focused checks for native score collection and index-based compilation."""

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scoring"))

from compile_scores import compile_scores
from ml_score_common import collect_method


def write_csv(path, rows, columns):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def write_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data), encoding="utf-8")


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


class NativeScoresTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.context = Path(self.temp.name) / "ctx_test"
        self.selected = self.context / "initial_screening" / "selected.csv"
        write_csv(self.selected, [{"SMILES": "CC"}, {"SMILES": "CC"}], ["SMILES"])

    def method_dir(self, name):
        path = self.context / "fine_screening" / name
        (path / "output").mkdir(parents=True, exist_ok=True)
        return path

    def test_nesso_duplicate_smiles_and_terminal_missing_score(self):
        path = self.method_dir("Nesso")
        write_json(path / "output/job_0/predictions/ctx_test_0/affinity.json",
                   {"affinity_pred_value": -1.2, "affinity_probability_binary": 0.8})
        (path / "output/job_ids.txt").write_text("1\n2\n", encoding="utf-8")
        (path / "output/job_results.tsv").write_text(
            "t\t1\tctx_test\tbatch_0\tsuccess\t1\t1\n"
            "t\t2\tctx_test\tbatch_1\tsuccess\t1\t1\n", encoding="utf-8")
        output, counts = collect_method("Nesso", path / "output", self.selected, path)
        rows = read_csv(output)
        self.assertEqual([r["status"] for r in rows], ["scored", "missing_score"])
        self.assertEqual([r["SMILES"] for r in rows], ["CC", "CC"])
        self.assertEqual(counts, {"scored": 1, "missing_score": 1})
        compiled, methods = compile_scores(self.context)
        self.assertEqual(methods, ["Nesso"])
        self.assertEqual(len(compiled), 2)
        self.assertEqual(compiled["nesso_status"].tolist(), ["scored", "missing_score"])

    def test_protenix_best_sample_and_tie_break(self):
        path = self.method_dir("Protenix")
        base = path / "output/job_0/ctx_test_0/seed_101/predictions"
        score = {"ranking_score": 0.9, "iptm": 0.7, "ptm": 0.8,
                 "chain_ptm": [0.8, 0.5], "chain_pair_iptm": [[0, 0.6], [0.6, 0]]}
        write_json(base / "ctx_test_0_summary_confidence_sample_1.json", score)
        write_json(base / "ctx_test_0_summary_confidence_sample_0.json", score)
        write_json(base / "ctx_test_0_summary_confidence_sample_2.json",
                   {**score, "ranking_score": 0.2})
        output, _ = collect_method("Protenix", path / "output", self.selected, path)
        rows = read_csv(output)
        self.assertEqual(rows[0]["sample"], "0")
        self.assertEqual(rows[0]["sample_count"], "3")
        self.assertEqual(rows[0]["ligand_chain_ptm"], "0.5")
        self.assertEqual(rows[1]["status"], "pending")

    def test_openfold_invalid_and_failed_job(self):
        path = self.method_dir("OpenFold3")
        base = path / "output/job_0/ctx_test_0/seed_42"
        write_json(base / "ctx_test_0_seed_42_sample_1_confidences_aggregated.json",
                   {"sample_ranking_score": "bad"})
        write_json(base / "ctx_test_0_seed_42_sample_2_confidences_aggregated.json",
                   {"sample_ranking_score": 0.8, "chain_ptm": {"A": 0.9, "Z": 0.5},
                    "chain_pair_iptm": {"(A, Z)": 0.6}})
        (path / "output/job_ids.txt").write_text("1\n2\n", encoding="utf-8")
        (path / "output/job_results.tsv").write_text(
            "t\t1\tctx_test\tbatch_0\tsuccess\t1\t1\n"
            "t\t2\tctx_test\tbatch_1\tfailed\t0\t1\n", encoding="utf-8")
        output, _ = collect_method("OpenFold3", path / "output", self.selected, path)
        rows = read_csv(output)
        self.assertEqual(rows[0]["status"], "scored")
        self.assertEqual(rows[0]["sample"], "2")
        self.assertEqual(rows[0]["sample_count"], "1")
        self.assertEqual(rows[0]["ligand_protein_pair_iptm_max"], "0.6")
        self.assertEqual(rows[1]["status"], "failed_job")

    def test_grouped_ligands_are_not_misassigned(self):
        path = self.method_dir("Nesso")
        manifest = path / "input/input_manifest.jsonl"
        manifest.parent.mkdir(parents=True, exist_ok=True)
        manifest.write_text(json.dumps({"compound_index": 0, "ligands": ["CC", "CC"]}) + "\n")
        output, _ = collect_method("Nesso", path / "output", self.selected, path)
        self.assertEqual([r["status"] for r in read_csv(output)],
                         ["unsupported_group", "unsupported_group"])

    def test_compiler_joins_existing_methods_by_row_index(self):
        fine = self.context / "fine_screening"
        write_csv(fine / "AF3/summary.csv", [
            {"folder": "job_1", "SMILES": "CC", "ranking_score": 0.7},
        ], ["folder", "SMILES", "ranking_score"])
        write_csv(fine / "Vina/input/input_part_0.csv", [
            {"ligand_description": "CC"}, {"ligand_description": "CC"},
        ], ["ligand_description"])
        write_csv(fine / "Vina/results.csv", [
            {"folder": 0, "idx": 0, "box": 0, "affinity": -6.0, "SMILES": "CC"},
            {"folder": 0, "idx": 0, "box": 1, "affinity": -8.0, "SMILES": "CC"},
            {"folder": 0, "idx": 1, "box": 0, "affinity": -7.0, "SMILES": "CC"},
        ], ["folder", "idx", "box", "affinity", "SMILES"])
        frame, methods = compile_scores(self.context)
        self.assertEqual(methods, ["AF3", "Vina"])
        self.assertEqual(len(frame), 2)
        self.assertTrue(str(frame.iloc[0]["af3_ranking_score"]) == "nan")
        self.assertEqual(frame["vina_affinity"].tolist(), [-8.0, -7.0])
        self.assertEqual(frame.iloc[1]["af3_ranking_score"], 0.7)

    def test_compiler_accepts_empty_vina_results_without_input_partitions(self):
        fine = self.context / "fine_screening"
        write_csv(fine / "Vina/results.csv", [],
                  ["folder", "idx", "box", "affinity", "SMILES"])
        frame, methods = compile_scores(self.context)
        self.assertEqual(methods, ["Vina"])
        self.assertEqual(len(frame), 2)
        self.assertTrue(frame["vina_affinity"].isna().all())

    def test_run_wide_merge_preserves_duplicate_smiles(self):
        run_root = self.context.parent
        write_csv(run_root / "predictions.csv", [
            {"prediction_id": "ctx_test__0", "context": "ctx_test", "SMILES": "CC"},
            {"prediction_id": "ctx_test__1", "context": "ctx_test", "SMILES": "CC"},
        ], ["prediction_id", "context", "SMILES"])
        script = Path(__file__).resolve().parents[1] / "scoring" / "collect_ml_run_scores.py"
        subprocess.run([sys.executable, str(script), "--run-root", str(run_root)],
                       check=True, capture_output=True, text=True)
        rows = read_csv(run_root / "compiled_scores.csv")
        self.assertEqual(len(rows), 2)
        self.assertEqual([row["prediction_id"] for row in rows],
                         ["ctx_test__0", "ctx_test__1"])
        self.assertTrue(all(row["nesso_status"] == "pending" for row in rows))


if __name__ == "__main__":
    unittest.main()
