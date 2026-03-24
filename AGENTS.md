# Repository Guidelines

## Project Structure & Module Organization
Core automation lives in `scripts/`. Use `run_full_pipeline.sh` for the 6-stage pipeline and `run_all_workflows.sh` for fine-screening orchestration only. Method job submitters use `run_*_batch.sh`; score collectors use `*_scores.py`.

`scripts/initial_screening/wrappers/` contains adapters for external tools (for example ConPLex and DrugLAMP). `scripts/md_pbsa/pbsa/` contains MD+PBSA helpers and MDP templates. Model implementations are vendored in `algos/` (`GraphDTA/`, `HMSA-DTI/`, `coldDTA/`) with method-specific dependency files.

Project documentation is in root `*.md` files. Runtime data is expected under `${TASK_ROOT}` on shared storage and should not be committed.

## Pipeline Stages
`run_full_pipeline.sh` executes these stages:

1. `make_input`: split sequence/smiles inputs into per-protein chunks.
2. `initial_screening`: submit GraphDTA, HMSA, ColdDTA, DrugLAMP, and ConPLex jobs.
3. `compile_summary`: aggregate method outputs and select top compounds.
4. `prepare_fine`: generate prefold/fine-screening inputs.
5. `fine_screening`: run AF3, Boltz2, RoseTTAFold, Vina, DiffDock, and MD+PBSA.
6. `collect_results`: parse outputs into per-method summary files.

Use `--start-from N` and `--stop-after N` to resume or limit stage execution.

## Software & Environment
Expected execution environment is HPC (AWS ParallelCluster-style) with SLURM (`sbatch`, `squeue`) and shared filesystem paths for `${TASK_ROOT}`.

Primary runtime stack:

- Bash orchestration scripts (`#!/bin/bash`, `set -e`).
- Conda-based Python environments (for example `conda activate HMSA` in input prep).
- Python tools used by scripts include `pandas`, `polars`, `rdkit`, `prody`, `openmm`, `pdbfixer`, and `vina`.
- Workflow methods include GraphDTA, HMSA-DTI, ColdDTA, AF3, Boltz2, RoseTTAFold, Vina, DiffDock, and MD+PBSA.

## Build, Test, and Development Commands
There is no single build step; this repository is script-driven.

- `./scripts/run_full_pipeline.sh --dry-run all`: validate stage wiring without submitting jobs.
- `./scripts/run_full_pipeline.sh --task-root /shared/task_1 --proteins "JAK1JH1 EGFR"`: run pipeline.
- `./scripts/run_full_pipeline.sh status`: stage-level status.
- `./scripts/run_all_workflows.sh status`: fine-screening workflow status.
- `./scripts/check_screening_progress.sh --task-root /shared/task_1 --proteins "JAK1JH1"`: method-level completion/failure summary.

## Coding Style & Naming Conventions
Follow existing conventions in `scripts/`:

- Bash: uppercase config vars (`TASK_ROOT`, `PROTEINS`), lowercase function names.
- Python: 4-space indentation and `snake_case` names.
- New script names should match existing patterns (`run_*.sh`, `*_scores.py`, `*_predict.py`).
- Prefer environment overrides (`MASTER_TASK_ROOT`, `MASTER_PROTEINS`) over new hardcoded paths.

## Testing Guidelines
Formal unit tests are not currently checked in. Minimum validation for changes:

1. `./scripts/run_full_pipeline.sh --dry-run all`
2. `./scripts/run_full_pipeline.sh status`
3. For parser/output changes, verify tokens and summary files under one real `${TASK_ROOT}/${PROTEIN}` run.

## Commit & Pull Request Guidelines
Keep commit subjects short and imperative (for example `fix run_full_pipeline.sh`, `add missing files`). In PRs, include changed stages/methods, exact validation commands run, and any SLURM/path assumptions.
