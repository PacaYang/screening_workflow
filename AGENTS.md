# Repository Guidelines

This is a branch designed specifically for GCP. Do not perform git add, git commit, git pull, git revert;

## Project Structure & Module Organization
- `scripts/pipeline/` contains the main orchestration layer; `run_full_pipeline.sh` is the primary entry point.
- `scripts/pipeline/stages/` holds stage-specific drivers (`stage1_make_input.sh` through `stage6_collect_results.sh`).
- `scripts/pipeline/lib/` provides shared Bash helpers for logging, state tracking, config loading, and SLURM job monitoring.
- Workflow logic is grouped by domain: `scripts/input/`, `scripts/initial_screening/`, `scripts/structure_prediction/`, `scripts/docking/`, `scripts/md_pbsa/`, `scripts/scoring/`, `scripts/monitoring/`, and `scripts/testing/`.

## Screening Stages & Software Used
- Stage 1 (`stage1_make_input.sh`) prepares split input tables from `Input/compounds_smiles.csv` and `Input/sequences.csv`.
- Stage 2 (`stage2_initial_screening.sh`) runs initial ML screening with GraphDTA, HMSA, ColdDTA, DrugLAMP, and ConPLex.
- Stage 3 (`stage3_compile_summary.sh`) merges model outputs and selects top compounds for fine screening.
- Stage 4 (`stage4_prepare_fine.sh`) writes method-specific inputs (AF3 JSON, Boltz2 YAML, and docking partitions).
- Stage 5 (`stage5_fine_screening.sh`) submits SLURM jobs for AF3, Boltz2, RoseTTAFold, Vina, DiffDock, and MD/PBSA.
- Stage 6 (`stage6_collect_results.sh`) normalizes outputs; default collection is streaming via `stage6_collect_results_streaming.sh`.
- Toolchain: Bash + Python (`pandas`, `numpy`, `argparse`), SLURM (`sbatch`, `squeue`), and external tools (AutoDock Vina, DiffDock, GROMACS, `gmx_MMPBSA`).

## Build, Test, and Development Commands
- Full pipeline:
  `bash scripts/pipeline/run_full_pipeline.sh --task-root /path/to/task --proteins "P1 P2" run`
- Status/resume:
  `bash scripts/pipeline/run_full_pipeline.sh --task-root /path/to/task status`
  `bash scripts/pipeline/run_full_pipeline.sh --task-root /path/to/task --resume-id <id> resume`
- Safe validation:
  `bash scripts/pipeline/run_full_pipeline.sh --task-root /path/to/task --dry-run run`
- Monitoring/smoke tests:
  `bash scripts/monitoring/check_screening_progress.sh`
  `bash scripts/testing/test_full_pipeline.sh`
  `sbatch scripts/testing/test_rosettafold_single.sh` (SLURM + GPU)

## Coding Style & Naming Conventions
- Bash conventions: use `#!/bin/bash`, `set -e`, uppercase config variables, and lowercase function names.
- Python conventions: 4-space indentation, `snake_case` identifiers, `argparse` for CLI arguments, and explicit output directories.
- Naming patterns: runners use `run_<workflow>.sh`, stages use `stage<index>_<action>.sh`, and scoring parsers use `<method>_scores.py`.

## Testing Guidelines
- This repo currently relies on integration/smoke testing rather than a unit-test framework or fixed coverage gate.
- For pipeline edits, run `--dry-run` first, then execute only affected stages with `--start-from` and `--stop-after`.
- Verify artifacts in method-specific output folders and state files (`.pipeline_state/`, `.collection_state/`).

