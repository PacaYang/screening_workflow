# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Repo Does

This is a multi-stage drug screening pipeline that runs on an HPC cluster (AWS ParallelCluster) via SLURM. It screens compound libraries against protein targets using a two-phase approach:

1. **Initial screening** — fast ML-based binding affinity prediction (GraphDTA, HMSA, ColdDTA, DrugLAMP, ConPLex)
2. **Fine screening** — physics-based and structure-prediction methods (Boltz2, AF3, RoseTTAFold, Vina docking, DiffDock, MD+PBSA)

All data lives on a shared FSx filesystem under a `TASK_ROOT` directory, organized per protein target.

## Running the Pipeline

The main entry point is `scripts/run_full_pipeline.sh`. It orchestrates 6 stages, submits SLURM jobs, and polls for completion.

```bash
# Full end-to-end run
./scripts/run_full_pipeline.sh --task-root /shared/task_1 --proteins "JAK1JH1 EGFR"

# Dry run to preview what would execute
./scripts/run_full_pipeline.sh --dry-run all

# Check progress across all stages
./scripts/run_full_pipeline.sh status

# Resume from a specific stage (e.g., after a failure)
./scripts/run_full_pipeline.sh --start-from 4 --task-root /shared/task_1 --proteins "JAK1JH1"

# Run only initial screening (stages 1-3)
./scripts/run_full_pipeline.sh --stop-after 3

# Skip specific methods
./scripts/run_full_pipeline.sh --skip-boltz2 --skip-md-pbsa
```

### Pipeline Stages

| Stage | Name | Description |
|-------|------|-------------|
| 1 | make_input | Split compound CSV into per-protein input chunks |
| 2 | initial_screening | Submit SLURM jobs for all 5 ML methods, wait for completion |
| 3 | compile_summary | Run `init_select_top.py` to select top compounds → `selected.csv` |
| 4 | prepare_fine | Write prefold inputs, run prefold SLURM jobs, write fine screening inputs |
| 5 | fine_screening | Submit Boltz2/AF3/Vina/DiffDock/MD+PBSA SLURM jobs; MD+PBSA waits on DiffDock |
| 6 | collect_results | Parse scores from outputs into per-method `summary.csv` files |

### Individual Workflow Scripts

Each fine-screening method has its own batch submission script:
- `scripts/run_boltz2_batch.sh` — submits 40 batch jobs
- `scripts/run_vina_batch.sh` — submits per-CSV-part jobs
- `scripts/run_diffdock_batch.sh` — submits per-CSV-part jobs
- `scripts/run_md_pbsa_batch.sh` — submits per-compound MD+PBSA jobs

Monitor progress:
```bash
squeue -u $USER
./scripts/monitor_boltz2_jobs.sh
./scripts/monitor_vina_jobs.sh
```

## Architecture

### Configuration Passing

Sub-scripts read configuration from environment variables set by the master script:
- `MASTER_TASK_ROOT` — task root directory
- `MASTER_PROTEINS` — space-separated protein list

Each sub-script falls back to hardcoded defaults if these are not set.

### Key Python Scripts

- `scripts/init_select_top.py` — selects top compounds from initial screening using Polars; seeds with HMSA ≥ threshold, fills remaining slots by best score across other methods
- `scripts/split_csv.py` — splits compound CSV into chunks for parallel SLURM jobs
- `scripts/boltz2_scores.py`, `af3_scores.py`, `vina_scores.py` — parse method-specific output formats into `summary.csv`
- `scripts/docking.py` — Vina docking logic
- `scripts/wrappers/run_conplex.py`, `run_druglamp.py`, `run_evidti.py` — wrappers for external tools

### Algorithm Implementations (`algos/`)

- `algos/GraphDTA/` — GNN-based DTA prediction; entry point `predict.py`, models in `models/`
- `algos/HMSA-DTI/` — HMSA binding prediction; entry point `predict_binding.py`
- `algos/coldDTA/` — ColdDTA prediction; entry point `predict.py`

Pre-trained model weights are stored alongside each algo (`model_para/` or `model/`).

### Completion Tokens

Stages use sentinel files to track completion and make re-runs idempotent:
- `finish.token` — stage 1 input preparation done
- `*.done` files — individual SLURM job completion markers
- `submitted_jobs.txt` — SLURM job IDs written by batch scripts, read by `wait_for_slurm_jobs()`

### Data Flow

```
Input/compounds_smiles.csv + sequences.csv
  → split into chunks → initial_screening/{GraphDTA,HMSA,ColdDTA,...}/prediction_*.csv
  → init_select_top.py → initial_screening/selected.csv
  → fine_screening/{Boltz2,AF3,Vina,PBSA}/...
  → score collection → fine_screening/{method}/summary.csv
```

## Branch Notes

Current branch is `AWS_parallelCluster`. The `main` branch is the base for PRs. When change the codes, do not perform any git operations such as git pull, git add, git commit, git revert, etc.
