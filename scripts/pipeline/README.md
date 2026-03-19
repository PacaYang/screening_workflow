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
