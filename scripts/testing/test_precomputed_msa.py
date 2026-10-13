#!/usr/bin/env python3
"""Offline smoke tests for the shared precomputed-MSA workflow."""

import gzip
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
MSA_HELPER = REPO_ROOT / "scripts/structure_prediction/colabfold_msa.py"
ML_GENERATOR = REPO_ROOT / "scripts/input/gen_ml_method_input.py"
BOLTZ_GENERATOR = REPO_ROOT / "scripts/input/gen_boltz_yaml.py"
BOLTZ_PYTHON = Path("/home/yangl_pacagen_com/miniconda3/envs/boltz_test/bin/python")
AF3_PYTHON = Path("/home/yangl_pacagen_com/miniconda3/envs/af3_test/bin/python")
OPENFOLD3_PYTHON = Path("/home/yangl_pacagen_com/miniconda3/envs/openfold3/bin/python")
PROTENIX_PYTHON = Path("/home/yangl_pacagen_com/miniconda3/envs/protenix/bin/python")
MSA_LAUNCHER = REPO_ROOT / "scripts/structure_prediction/run_colabfold_msa.sh"


class PrecomputedMsaWorkflowTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.sequences = self.root / "sequences.csv"
        self.selected = self.root / "selected.csv"
        self.sequences.write_text(
            "protein_name,sequence\nP1,\"ACDE,ACDE,MNPQ\"\n", encoding="utf-8"
        )
        self.selected.write_text("SMILES\nCCO\nCCN\n", encoding="utf-8")
        self.work = self.root / "work"
        self.msa_dir = self.root / "MSA"

    def tearDown(self):
        self.temporary.cleanup()

    def run_python(self, script, *arguments, check=True, executable=None):
        return subprocess.run(
            [str(executable or sys.executable), str(script), *map(str, arguments)],
            check=check,
            capture_output=True,
            text=True,
        )

    def prepare_manifest(self):
        self.run_python(
            MSA_HELPER,
            "prepare",
            "--sequences-csv",
            self.sequences,
            "--protein",
            "P1",
            "--work-dir",
            self.work,
        )
        plan = json.loads((self.work / "search_plan.json").read_text(encoding="utf-8"))
        self.assertEqual(len(plan["searches"]), 2)
        sequence_by_chain = {chain["chain_id"]: chain["sequence"] for chain in plan["chains"]}
        for search in plan["searches"]:
            output_dir = Path(search["output_dir"])
            output_dir.mkdir(parents=True)
            sequence = sequence_by_chain[search["representative_chain"]]
            (output_dir / "result.a3m").write_text(
                f"#{len(sequence)}\t1\n>query\n{sequence}\n>hit\n{sequence}\n",
                encoding="utf-8",
            )
        self.run_python(
            MSA_HELPER,
            "finalize",
            "--work-dir",
            self.work,
            "--msa-dir",
            self.msa_dir,
            "--colabfold-command",
            "fake_colabfold_batch --msa-only",
        )
        self.run_python(
            MSA_HELPER,
            "check",
            "--sequences-csv",
            self.sequences,
            "--protein",
            "P1",
            "--msa-dir",
            self.msa_dir,
        )
        return self.msa_dir / "manifest.json"

    def test_manifest_deduplicates_search_and_normalizes_each_chain(self):
        manifest_path = self.prepare_manifest()
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        self.assertEqual(manifest["format_version"], 2)
        self.assertEqual([entry["chain_id"] for entry in manifest["chains"]], ["A", "B", "C"])
        self.assertEqual(manifest["chains"][0]["representative_chain"], "A")
        self.assertEqual(manifest["chains"][1]["representative_chain"], "A")
        self.assertTrue(all(Path(entry["msa_path"]).is_file() for entry in manifest["chains"]))
        self.assertEqual(
            len({Path(entry["msa_path"]).parent for entry in manifest["chains"]}),
            2,
        )
        self.assertEqual(manifest["chains"][0]["msa_path"], manifest["chains"][1]["msa_path"])

    def test_method_inputs_reference_shared_a3m_files(self):
        manifest_path = self.prepare_manifest()
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        msa_paths = [entry["msa_path"] for entry in manifest["chains"]]

        for method, suffix in (
            ("AF3", ".json"),
            ("Protenix", ".protenix.json"),
            ("OpenFold3", ".openfold.json"),
        ):
            output_dir = self.root / method
            self.run_python(
                ML_GENERATOR,
                "--method",
                method,
                "--mode",
                "complex",
                "--output-dir",
                output_dir,
                "--input-csv",
                self.sequences,
                "--protein-name",
                "P1",
                "--smiles-file",
                self.selected,
                "--msa-manifest",
                manifest_path,
            )
            payload = json.loads((output_dir / f"P1_0{suffix}").read_text(encoding="utf-8"))
            if method == "AF3":
                self.assertNotIn("screening", payload)
                proteins = [entry["protein"] for entry in payload["sequences"] if "protein" in entry]
                self.assertEqual([entry["unpairedMsaPath"] for entry in proteins], msa_paths)
                self.assertTrue(all(entry["pairedMsa"] == "" for entry in proteins))
                subprocess.run(
                    [
                        str(AF3_PYTHON),
                        "-c",
                        "import sys; from pathlib import Path; "
                        "sys.path.insert(0, '/home/yangl_pacagen_com/Applications/alphafold3/src'); "
                        "from alphafold3.common import folding_input; "
                        "folding_input.load_fold_inputs_from_path(Path(sys.argv[1]))",
                        str(output_dir / "P1_0.json"),
                    ],
                    check=True,
                )
            elif method == "Protenix":
                proteins = [entry["proteinChain"] for entry in payload[0]["sequences"] if "proteinChain" in entry]
                self.assertEqual([entry["unpairedMsaPath"] for entry in proteins], msa_paths)
                subprocess.run(
                    [
                        str(PROTENIX_PYTHON),
                        "-c",
                        "import json, sys; from pathlib import Path; "
                        "from protenix.data.constants import PROTEIN_CHAIN; "
                        "from protenix.data.msa.msa_utils import RawMsa; "
                        "payload = json.loads(Path(sys.argv[1]).read_text())[0]; "
                        "chains = [x['proteinChain'] for x in payload['sequences'] "
                        "if 'proteinChain' in x]; "
                        "[RawMsa.from_a3m(c['sequence'], PROTEIN_CHAIN, "
                        "Path(c['unpairedMsaPath']).read_text()) for c in chains]",
                        str(output_dir / "P1_0.protenix.json"),
                    ],
                    check=True,
                )
            else:
                query = payload["queries"]["P1_0"]
                proteins = [chain for chain in query["chains"] if chain["molecule_type"] == "protein"]
                self.assertTrue(query["use_msas"])
                self.assertTrue(query["use_main_msas"])
                self.assertFalse(query["use_paired_msas"])
                self.assertEqual([entry["main_msa_file_paths"][0] for entry in proteins], msa_paths)
                subprocess.run(
                    [
                        str(OPENFOLD3_PYTHON),
                        "-c",
                        "import sys; from pathlib import Path; "
                        "from openfold3.projects.of3_all_atom.config.inference_query_format "
                        "import InferenceQuerySet; "
                        "from openfold3.core.config.msa_pipeline_configs import "
                        "MsaSampleProcessorInputInference; "
                        "from openfold3.core.data.io.sequence.msa import "
                        "MsaSampleParserInference; "
                        "from openfold3.projects.of3_all_atom.config.dataset_config_components "
                        "import MSASettings; "
                        "iqs = InferenceQuerySet.from_json(Path(sys.argv[1])); "
                        "query = next(iter(iqs.queries.values())); "
                        "msa_input = MsaSampleProcessorInputInference.create_from_inference_query_entry(query); "
                        "parser = MsaSampleParserInference(MSASettings()); "
                        "maps = parser.create_maps(msa_input); "
                        "assert len(set(maps.chain_id_to_rep_id.values())) == 2; "
                        "assert maps.chain_id_to_rep_id['A'] == maps.chain_id_to_rep_id['B']; "
                        "assert maps.chain_id_to_rep_id['A'] != maps.chain_id_to_rep_id['C']; "
                        "maps = parser.parse_msas(maps); "
                        "assert len(maps.rep_id_to_main_msa) == 2",
                        str(output_dir / "P1_0.openfold.json"),
                    ],
                    check=True,
                )

        boltz_dir = self.root / "Boltz2"
        self.run_python(
            BOLTZ_GENERATOR,
            "--output",
            boltz_dir,
            "--protein-file",
            self.sequences,
            "--protein-name",
            "P1",
            "--name-col",
            "protein_name",
            "--smiles-path",
            self.selected,
            "--msa-manifest",
            manifest_path,
            executable=BOLTZ_PYTHON,
        )
        boltz_text = (boltz_dir / "0.yaml").read_text(encoding="utf-8")
        for msa_path in msa_paths:
            self.assertIn(f"msa: {msa_path}", boltz_text)
        subprocess.run(
            [
                str(BOLTZ_PYTHON),
                "-c",
                "import json, sys; from pathlib import Path; "
                "from rdkit.Chem import AllChem, Descriptors; "
                "AllChem.Descriptors = Descriptors; "
                "from boltz.data.parse.a3m import parse_a3m; "
                "from boltz.data.parse.yaml import parse_yaml; "
                "yaml_path = Path(sys.argv[1]); "
                "parse_yaml(yaml_path, {}, Path('/home/yangl_pacagen_com/.boltz/mols'), boltz2=True); "
                "manifest = json.loads(Path(sys.argv[2]).read_text()); "
                "[parse_a3m(Path(entry['msa_path']), taxonomy=None, max_seqs=16384) "
                "for entry in manifest['chains']]",
                str(boltz_dir / "0.yaml"),
                str(manifest_path),
            ],
            check=True,
        )

    def test_changed_sequence_invalidates_cache(self):
        self.prepare_manifest()
        self.sequences.write_text(
            "protein_name,sequence\nP1,\"ACDF,ACDE,MNPQ\"\n", encoding="utf-8"
        )
        result = self.run_python(
            MSA_HELPER,
            "check",
            "--sequences-csv",
            self.sequences,
            "--protein",
            "P1",
            "--msa-dir",
            self.msa_dir,
            check=False,
        )
        self.assertNotEqual(result.returncode, 0)

    def test_gzip_a3m_output_is_normalized(self):
        self.run_python(
            MSA_HELPER,
            "prepare",
            "--sequences-csv",
            self.sequences,
            "--protein",
            "P1",
            "--work-dir",
            self.work,
        )
        plan = json.loads((self.work / "search_plan.json").read_text(encoding="utf-8"))
        sequence_by_chain = {chain["chain_id"]: chain["sequence"] for chain in plan["chains"]}
        for search in plan["searches"]:
            output_dir = Path(search["output_dir"])
            output_dir.mkdir(parents=True)
            sequence = sequence_by_chain[search["representative_chain"]]
            with gzip.open(output_dir / "result.a3m.gz", "wt", encoding="utf-8") as handle:
                handle.write(f">query\n{sequence}\n>hit\n{sequence}\n")
        self.run_python(
            MSA_HELPER,
            "finalize",
            "--work-dir",
            self.work,
            "--msa-dir",
            self.msa_dir,
            "--colabfold-command",
            "fake_colabfold_batch --msa-only",
        )
        manifest = json.loads((self.msa_dir / "manifest.json").read_text(encoding="utf-8"))
        for entry in manifest["chains"]:
            msa_path = Path(entry["msa_path"])
            self.assertEqual(msa_path.suffix, ".a3m")
            self.assertTrue(msa_path.read_text(encoding="utf-8").startswith(">query\n"))

    def test_launcher_generates_valid_slurm_script(self):
        task_root = self.root / "task"
        (task_root / "Input").mkdir(parents=True)
        (task_root / "Input/sequences.csv").write_text(
            self.sequences.read_text(encoding="utf-8"), encoding="utf-8"
        )
        fake_colabfold = self.root / "fake_colabfold_batch.py"
        fake_colabfold.write_text(
            "#!/usr/bin/env python3\n"
            "import sys\n"
            "from pathlib import Path\n"
            "source, output = Path(sys.argv[1]), Path(sys.argv[2])\n"
            "sequence = ''.join(source.read_text().splitlines()[1:])\n"
            "output.mkdir(parents=True, exist_ok=True)\n"
            "(output / 'result.a3m').write_text(f'>query\\n{sequence}\\n>hit\\n{sequence}\\n')\n",
            encoding="utf-8",
        )
        fake_colabfold.chmod(0o755)
        environment = os.environ.copy()
        environment["MASTER_COLABFOLD_ENV"] = ""
        environment["MASTER_COLABFOLD_BATCH"] = str(fake_colabfold)
        subprocess.run(
            [
                "bash",
                str(MSA_LAUNCHER),
                "--task-root",
                str(task_root),
                "--proteins",
                "P1",
                "--dry-run",
            ],
            check=True,
            capture_output=True,
            text=True,
            env=environment,
        )
        job_script = task_root / "P1/fine_screening/MSA/logs/slurm_colabfold_msa.sh"
        self.assertTrue(job_script.is_file())
        subprocess.run(["bash", "-n", str(job_script)], check=True)
        environment["SLURM_JOB_ID"] = "offline-test"
        subprocess.run(["bash", str(job_script)], check=True, env=environment)
        self.run_python(
            MSA_HELPER,
            "check",
            "--sequences-csv",
            task_root / "Input/sequences.csv",
            "--protein",
            "P1",
            "--msa-dir",
            task_root / "P1/fine_screening/MSA",
        )


if __name__ == "__main__":
    unittest.main()
