# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.
This is the GCP branch of the repo, specifically used for GCP. 
When change the codes, do not perform any git operations, eg, git add, git commit, git revert, git pull, etc.

## What This Repo Does

Multi-stage computational drug screening pipeline. Screens compound libraries against protein targets using:
1. Fast ML-based binding affinity models (initial screening: ColdDTA, HMSA, GraphDTA, DrugLAMP, ConPLex)
2. Structure prediction + docking + MD/PBSA (fine screening: AF3, Boltz2, RoseTTAFold, Vina, DiffDock/PBSA)

Runs on a GCP cluster with SLURM job scheduling.

## Running the Pipeline

```bash
# Full pipeline (streaming collection is default)
bash scripts/pipeline/run_full_pipeline.sh \
  --task-root /path/to/task \
  --proteins "PROTEIN1 PROTEIN2" \
  run

# Resume from a specific stage
bash scripts/pipeline/run_full_pipeline.sh \
  --task-root /path/to/task \
  --proteins "PROTEIN1" \
  --start-from 5 \
  run

# Check status
bash scripts/pipeline/run_full_pipeline.sh --task-root /path/to/task status

# Dry run
bash scripts/pipeline/run_full_pipeline.sh --task-root /path/to/task --dry-run run
```

Key flags: `--skip-af3`, `--skip-boltz2`, `--skip-vina`, `--skip-diffdock`, `--skip-rosettafold`, `--skip-pbsa`, `--target-n N` (compounds for fine screening, default 3000), `--batch-collection` (wait for all jobs before collecting).

## Pipeline Stages

| Stage | Script | Purpose |
|-------|--------|---------|
| 1 | `stages/stage1_make_input.sh` | Generate input CSVs from SMILES + sequences |
| 2 | `stages/stage2_initial_screening.sh` | Run ML models (GraphDTA, HMSA, ColdDTA, DrugLAMP, ConPLex) |
| 3 | `stages/stage3_compile_summary.sh` | Aggregate initial screening results, select top N |
| 4 | `stages/stage4_prepare_fine.sh` | Prepare inputs for fine screening |
| 5 | `stages/stage5_fine_screening.sh` | AF3, Boltz2, RosettaFold, Vina, DiffDock, MD/PBSA (parallel SLURM jobs) |
| 6 | `stages/stage6_collect_results.sh` | Collect and compile results (streaming or batch) |

## Architecture

**Entry point:** `scripts/pipeline/run_full_pipeline.sh` — orchestrates all stages, delegates to `scripts/pipeline/stages/stage*.sh`.

**Shared libraries** in `scripts/pipeline/lib/`:
- `logger.sh` — logging
- `state_manager.sh` — JSON state in `.pipeline_state/` (enables resume)
- `job_monitor.sh` — SLURM job tracking
- `collection_state.sh` — streaming collection state per method

**Scoring scripts** in `scripts/scoring/` parse raw outputs from each method into CSVs. `compile_scores.py` merges all method results.

**Streaming collection** (default): results are collected incrementally as SLURM jobs finish. State tracked in `${TASK_ROOT}/${PROTEIN}/fine_screening/.collection_state/{AF3,Boltz2,Vina,PBSA}.json`. AF3 and Boltz2 support incremental collection; Vina and PBSA are pending.

## Expected Input Layout

```
${TASK_ROOT}/Input/
├── compounds_smiles.csv
├── sequences.csv
└── protein_file/${PROTEIN}/
    ├── ${PROTEIN}.pdb
    ├── ${PROTEIN}.gro
    └── system_EM.top
```

## Key Paths

- Model weights: `/home/yangl_pacagen_com/Applications/model_weights`

## Debugging Collection Issues

```bash
# Check collection state for a method
cat ${TASK_ROOT}/${PROTEIN}/fine_screening/.collection_state/AF3.json | jq .

# Reset and re-run collection
rm -rf ${TASK_ROOT}/${PROTEIN}/fine_screening/.collection_state/
bash scripts/pipeline/stages/stage6_collect_results_streaming.sh \
  --task-root ${TASK_ROOT} --proteins ${PROTEIN} --max-iterations 1
```
