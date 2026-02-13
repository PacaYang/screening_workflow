# Screening Workflow Automation Scripts - Summary

## All Scripts Created

### Master Control Script
- **`run_all_workflows.sh`** - Master orchestration script that controls all workflows

### Individual Workflow Scripts
1. **`run_boltz2_batch.sh`** - Boltz2 structure prediction (40 batches)
2. **`run_vina_batch.sh`** - AutoDock Vina docking (per CSV part)
3. **`run_diffdock_batch.sh`** - DiffDock docking (per CSV part)
4. **`run_md_pbsa_batch.sh`** - Combined MD + PBSA analysis (per compound)

### Monitoring Scripts
1. **`monitor_boltz2_jobs.sh`** - Monitor Boltz2 batch progress
2. **`monitor_vina_jobs.sh`** - Monitor Vina docking progress
3. **`monitor_diffdock_jobs.sh`** - Monitor DiffDock progress
4. **`monitor_md_pbsa_jobs.sh`** - Monitor MD+PBSA progress

### Documentation
1. **`BOLTZ2_BATCH_README.md`** - Boltz2 workflow documentation
2. **`MASTER_WORKFLOW_README.md`** - Master script documentation

## Quick Reference

### Run Everything
```bash
./run_all_workflows.sh
```

### Run Individual Workflows
```bash
./run_all_workflows.sh boltz2      # Boltz2 only
./run_all_workflows.sh vina        # Vina only
./run_all_workflows.sh diffdock    # DiffDock only
./run_all_workflows.sh md_pbsa     # MD+PBSA only
```

### Custom Configuration
```bash
./run_all_workflows.sh --task-root /my/data --proteins "JAK1 JAK2"
```

### Check Status
```bash
./run_all_workflows.sh status
```

### Monitor Progress
```bash
./monitor_boltz2_jobs.sh
./monitor_vina_jobs.sh
./monitor_diffdock_jobs.sh
./monitor_md_pbsa_jobs.sh
```

## Configuration Hierarchy

The scripts support configuration in three ways (priority order):

1. **Master script command line** (highest priority)
   ```bash
   ./run_all_workflows.sh --task-root /data --proteins "P1 P2"
   ```

2. **Environment variables**
   ```bash
   export MASTER_TASK_ROOT="/data"
   export MASTER_PROTEINS="P1 P2"
   ./run_boltz2_batch.sh
   ```

3. **Default values in scripts** (lowest priority)
   - Edit each script to change defaults
   - Located at top of each script

## Workflow Dependencies

```
Initial Screening (Snakemake)
    ↓
    ├─→ Boltz2 (independent)
    ├─→ Vina (independent)
    └─→ DiffDock
          ↓
          └─→ MD + PBSA
```

- **Boltz2** and **Vina** can run independently
- **DiffDock** requires Vina input preparation (split CSV files)
- **MD+PBSA** requires DiffDock output (rank1.sdf files)
  - Jobs skip gracefully if SDF files are missing

## Key Features

### Master Script (`run_all_workflows.sh`)
✓ Single command to run all workflows
✓ Selective workflow execution
✓ Centralized configuration
✓ Prerequisites checking
✓ Overall status reporting
✓ Help and usage information

### Boltz2 (`run_boltz2_batch.sh`)
✓ 40 batch jobs for parallel processing
✓ Automatic job range calculation
✓ Local /tmp storage for speed
✓ Result compression (tar.gz)
✓ EC2 instance type control (g5.xlarge)

### Vina (`run_vina_batch.sh`)
✓ Per-part job submission
✓ Automatic docking box extraction
✓ CPU-based docking (configurable GPU)
✓ Independent part processing

### DiffDock (`run_diffdock_batch.sh`)
✓ Per-part job submission
✓ GPU-accelerated docking
✓ Uses same input as Vina
✓ Flexible error handling

### MD+PBSA (`run_md_pbsa_batch.sh`)
✓ Combined MD and PBSA in single job
✓ Sequential execution (MD → PBSA)
✓ Graceful skipping of missing inputs
✓ Per-compound job submission
✓ Comprehensive logging

## SLURM Configuration

All scripts include:
- Partition specification
- Time limits
- Memory allocation
- CPU/GPU requirements
- **EC2 instance constraint** (`--constraint=g5.xlarge`)
- Job naming for easy identification
- Separate stdout/stderr logs

## File Structure

```
screening_workflow/
├── run_all_workflows.sh              # Master script
├── run_boltz2_batch.sh              # Boltz2 automation
├── run_vina_batch.sh                # Vina automation
├── run_diffdock_batch.sh            # DiffDock automation
├── run_md_pbsa_batch.sh             # MD+PBSA automation
├── monitor_boltz2_jobs.sh           # Boltz2 monitoring
├── monitor_vina_jobs.sh             # Vina monitoring
├── monitor_diffdock_jobs.sh         # DiffDock monitoring
├── monitor_md_pbsa_jobs.sh          # MD+PBSA monitoring
├── BOLTZ2_BATCH_README.md           # Boltz2 documentation
├── MASTER_WORKFLOW_README.md        # Master documentation
└── AUTOMATION_SUMMARY.md            # This file
```

## Common Tasks

### Submit all workflows for one protein
```bash
./run_all_workflows.sh --proteins "JAK1JH1"
```

### Submit all workflows for multiple proteins
```bash
./run_all_workflows.sh --proteins "JAK1 JAK2 JAK3"
```

### Run only structure prediction
```bash
./run_all_workflows.sh boltz2
```

### Run only docking
```bash
./run_all_workflows.sh --skip-boltz2 --skip-md-pbsa
```

### Run only binding affinity calculations
```bash
./run_all_workflows.sh md_pbsa
```

### Monitor all workflows
```bash
# Quick status check
./run_all_workflows.sh status

# Detailed monitoring
./monitor_boltz2_jobs.sh
./monitor_vina_jobs.sh
./monitor_diffdock_jobs.sh
./monitor_md_pbsa_jobs.sh

# SLURM queue
squeue -u $USER
```

### Cancel all jobs for a workflow
```bash
scancel -n boltz2_JAK1JH1_*
scancel -n vina_JAK1JH1_*
scancel -n diffdock_JAK1JH1_*
scancel -n md_pbsa_JAK1JH1_*
```

## Updates Made to Individual Scripts

All four workflow scripts now:
1. Check `MASTER_TASK_ROOT` environment variable (falls back to default)
2. Check `MASTER_PROTEINS` environment variable (falls back to default)
3. Can be called standalone or via master script
4. Maintain backward compatibility with direct usage

## Next Steps

1. **Configure paths**: Edit `run_all_workflows.sh` lines 12-16
2. **Make executable**: `chmod +x *.sh`
3. **Test**: Run `./run_all_workflows.sh help`
4. **Verify prerequisites**: Run `./run_all_workflows.sh status`
5. **Submit workflows**: Run `./run_all_workflows.sh`

## Support

For detailed documentation, see:
- `MASTER_WORKFLOW_README.md` - Master script usage
- `BOLTZ2_BATCH_README.md` - Boltz2 workflow details
- Individual script headers - Configuration options
