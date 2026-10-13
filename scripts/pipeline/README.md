# Modular Pipeline Architecture

This directory contains the modular implementation of the drug screening pipeline.

## Structure

```
scripts/pipeline/
├── run_full_pipeline.sh          # Master controller (NEW)
├── run_af3_first_pipeline.sh     # AF3-first controller (no initial screening)
├── run_full_pipeline.sh.backup   # Original monolithic script (backup)
├── lib/                          # Shared libraries
│   ├── logger.sh                 # Logging functions
│   ├── pipeline_utils.sh         # Common utilities
│   ├── state_manager.sh          # JSON state tracking
│   ├── job_monitor.sh            # SLURM job monitoring
│   └── config_loader.sh          # YAML configuration support
└── stages/                       # Stage scripts
    ├── stage1_make_input.sh
    ├── stage2_initial_screening.sh
    ├── stage3_compile_summary.sh
    ├── stage4_prepare_fine.sh
    ├── stage5_fine_screening.sh
    └── stage6_collect_results.sh
```

## Features

- **Modular Design**: Each stage is a separate script with clear inputs/outputs
- **State Management**: JSON-based progress tracking in `.pipeline_state/`
- **Resume Capability**: Resume interrupted pipelines from any stage
- **Enhanced Logging**: Structured logs in `.pipeline_logs/`
- **Parallel Execution**: Independent jobs run simultaneously
- **Backward Compatible**: Same CLI interface as original script

## Usage

### Run full pipeline
```bash
./run_full_pipeline.sh --task-root /data/screen --proteins "JAK1JH1 EGFR" run
```

### Run AF3-first pipeline (no initial screening)
```bash
./run_af3_first_pipeline.sh --task-root /data/screen --proteins "JAK1JH1 EGFR" run
```

### Check status
```bash
./run_full_pipeline.sh --task-root /data/screen status
```

### Resume interrupted pipeline
```bash
./run_full_pipeline.sh --task-root /data/screen resume
```

### Dry run
```bash
./run_full_pipeline.sh --task-root /data/screen --dry-run run
```

### Reuse one ColabFold MSA across structure methods

```bash
./run_full_pipeline.sh \
  --task-root /data/screen \
  --proteins "JAK1JH1 EGFR" \
  --precompute-msa \
  run
```

This opt-in mode searches the public ColabFold server once per unique protein-
chain sequence, then supplies the resulting A3M to AF3, Boltz2, Protenix, and
OpenFold3.
The normalized alignments and a sequence-hash manifest are stored in
`${TASK_ROOT}/${PROTEIN}/fine_screening/MSA`. Identical chains share one search
and one stable A3M path; different sequences receive distinct paths.
`colabfold_batch` must be available in the `colabfold` conda environment, or
configured with
`MASTER_COLABFOLD_ENV` and `MASTER_COLABFOLD_BATCH`.

Because the default mode uses the public ColabFold server, protein sequences
leave the cluster and server availability/rate limits apply. The option cannot
be combined with `--template-mode`. RoseTTAFold and Nesso keep their existing
prefold behavior.

The MSA search can also be submitted independently:

```bash
MASTER_TASK_ROOT=/data/screen MASTER_PROTEINS="JAK1JH1 EGFR" \
  bash scripts/structure_prediction/run_colabfold_msa.sh
```

### Run specific stages
```bash
./run_full_pipeline.sh --start-from 4 --stop-after 5 run
```

## State Files

Pipeline state is tracked in JSON format:
- Location: `${TASK_ROOT}/.pipeline_state/${pipeline_id}.json`
- Contains: stage status, job IDs, timestamps, configuration
- Used for: resume capability, progress tracking, debugging

## Migration from Old Script

The original script is backed up as `run_full_pipeline.sh.backup`. The new modular version maintains the same CLI interface, so existing workflows should work without changes.

## Benefits

1. **Maintainability**: ~200 lines per component vs 850 lines monolithic
2. **Testability**: Each stage can be tested independently
3. **Debuggability**: Easier to isolate and fix issues
4. **Extensibility**: Simple to add new stages or methods
5. **Observability**: JSON state + structured logging
6. **Reliability**: Robust error handling and resume capability

## AF3-like ML methods

The fine-screening stages also support NESSO, Protenix, and OpenFold3. They
use the same protein/SMILES inputs as AF3 and emit the native input expected by
each installed package (`.yaml`, `.protenix.json`, and `.openfold.json`). Native
outputs and a per-batch ledger are written under `fine_screening/Nesso`,
`fine_screening/Protenix`, and `fine_screening/OpenFold3`. Stage 6 extracts
native scores into each method's `summary.csv` and compiles available method
summaries into `fine_screening/compiled_scores.csv`. See
`scripts/scoring/ML_METHOD_SCORING.md` for score fields and standalone usage.

```bash
bash run_full_pipeline.sh \
  --task-root /data/screening \
  --proteins "JAK1JH1 EGFR" \
  --skip-af3 --skip-boltz2 --skip-rosettafold \
  run
```

The package environments default to `nesso`, `protenix`, and `openfold3` under
the shared Miniconda installation. Persistent checkpoints and generated caches
default to `/home/yangl_pacagen_com/Applications/model_weights/{Nesso,Protenix,OpenFold3}`.
These values can be overridden with `MASTER_NESSO_*`, `MASTER_PROTENIX_*`, and
`MASTER_OPENFOLD3_*` variables (`ROOT`, `ENV`, `EXE`, `CMD`, `WEIGHTS`, and
`DB`).
The command variables accept `{root}`, `{input}`, `{output}`, `{method}`,
`{exe}`, `{weights}`, and `{database}` placeholders.
