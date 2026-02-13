# Master Workflow Automation Script

This master script orchestrates all four screening workflows: Boltz2, Vina, DiffDock, and MD+PBSA.

## Overview

The `run_all_workflows.sh` script provides a unified interface to:
- Run all workflows with a single command
- Run individual workflows selectively
- Configure task root and protein list from command line
- Check overall status of all workflows
- Pass configuration to all sub-scripts automatically

## Quick Start

```bash
# Run all workflows with default settings
./run_all_workflows.sh

# Run with custom task root and proteins
./run_all_workflows.sh --task-root /path/to/data --proteins "PROTEIN1 PROTEIN2 PROTEIN3"

# Run only specific workflows
./run_all_workflows.sh boltz2
./run_all_workflows.sh vina
./run_all_workflows.sh diffdock
./run_all_workflows.sh md_pbsa

# Check status of all workflows
./run_all_workflows.sh status

# Show help
./run_all_workflows.sh help
```

## Configuration

Edit the script to set default values (lines 12-16):

```bash
# Task root directory
TASK_ROOT="/home/ubuntu/snake_test"

# Protein list (space-separated)
PROTEINS="JAK1JH1"
```

Or override from command line:

```bash
./run_all_workflows.sh --task-root /my/data --proteins "JAK1 JAK2 JAK3"
```

## Commands

### Run All Workflows (default)

```bash
./run_all_workflows.sh
# or explicitly
./run_all_workflows.sh all
```

Runs all enabled workflows in sequence:
1. Boltz2 (40 batches)
2. Vina (per-part jobs)
3. DiffDock (per-part jobs)
4. MD+PBSA (per-compound jobs)

### Run Individual Workflows

```bash
# Run only Boltz2
./run_all_workflows.sh boltz2

# Run only Vina
./run_all_workflows.sh vina

# Run only DiffDock
./run_all_workflows.sh diffdock

# Run only MD+PBSA
./run_all_workflows.sh md_pbsa
```

### Skip Workflows

```bash
# Run all except DiffDock and MD+PBSA
./run_all_workflows.sh --skip-diffdock --skip-md-pbsa

# Run only Boltz2 and Vina
./run_all_workflows.sh --skip-diffdock --skip-md-pbsa
```

### Check Status

```bash
./run_all_workflows.sh status
```

Shows:
- Completion status for each workflow
- Number of completed jobs vs total
- SLURM queue status (running and pending jobs)

Example output:
```
Workflow Status
==========================================

Protein: JAK1JH1
----------------------------------------
  Boltz2: 35/40 batches completed
  Vina: 8/10 parts completed
  DiffDock: 6/10 parts completed
  MD: 45/102 completed
  PBSA: 40/102 completed

SLURM Queue Status:
----------------------------------------
Running jobs:
  boltz2_JAK1JH1_b36 RUNNING
  vina_JAK1JH1_part_8 RUNNING
  md_pbsa_JAK1JH1_50 RUNNING

Pending jobs:
  boltz2_JAK1JH1_b37 PENDING
  boltz2_JAK1JH1_b38 PENDING
```

## Options

### --task-root DIR
Set the task root directory. All protein subdirectories are expected under this path.

```bash
./run_all_workflows.sh --task-root /mnt/data/screening
```

### --proteins LIST
Set the protein list (space-separated). The master script will process all listed proteins.

```bash
./run_all_workflows.sh --proteins "JAK1JH1 JAK2 JAK3"
```

### --skip-[workflow]
Skip specific workflows:
- `--skip-boltz2` - Skip Boltz2
- `--skip-vina` - Skip Vina
- `--skip-diffdock` - Skip DiffDock
- `--skip-md-pbsa` - Skip MD+PBSA

```bash
# Run only Vina and DiffDock
./run_all_workflows.sh --skip-boltz2 --skip-md-pbsa
```

## How It Works

### Configuration Passing

The master script exports configuration as environment variables:
- `MASTER_TASK_ROOT` - Task root directory
- `MASTER_PROTEINS` - Space-separated protein list

Each sub-script (`run_boltz2_batch.sh`, `run_vina_batch.sh`, etc.) checks for these variables:

```bash
# In sub-scripts
TASK_ROOT="${MASTER_TASK_ROOT:-/home/ubuntu/snake_test}"  # Use env var or default

get_proteins() {
    if [ -n "$MASTER_PROTEINS" ]; then
        echo "$MASTER_PROTEINS"  # Use env var
    else
        echo "JAK1JH1"  # Use default
    fi
}
```

### Prerequisites Check

Before running workflows, the script checks:
- Task root directory exists
- Protein directories exist
- Initial screening results (`selected.csv`) exist

### Workflow Execution

Each workflow is executed by calling its respective script:
1. `run_boltz2_batch.sh` - Submits 40 batch jobs
2. `run_vina_batch.sh` - Submits jobs for each CSV part
3. `run_diffdock_batch.sh` - Submits jobs for each CSV part
4. `run_md_pbsa_batch.sh` - Submits combined MD+PBSA jobs for each compound

## Directory Structure

Expected directory structure:

```
${TASK_ROOT}/
├── Input/
│   ├── sequences.csv
│   └── protein_file/
│       └── ${PROTEIN}/
│           ├── ${PROTEIN}.pdb
│           ├── ${PROTEIN}.gro
│           └── system_EM.top
└── ${PROTEIN}/
    ├── initial_screening/
    │   └── selected.csv
    └── fine_screening/
        ├── Boltz2/
        ├── Vina/
        └── PBSA/
            ├── DiffDock/
            └── PBSA/
                ├── MD/
                └── PBSA/
```

## Monitoring

After submission, use individual monitoring scripts:

```bash
# Boltz2 progress
./monitor_boltz2_jobs.sh

# Vina progress
./monitor_vina_jobs.sh

# DiffDock progress
./monitor_diffdock_jobs.sh

# MD+PBSA progress
./monitor_md_pbsa_jobs.sh

# Or check overall status
./run_all_workflows.sh status

# Or use SLURM commands
squeue -u $USER
squeue -u $USER | grep boltz2
squeue -u $USER | grep vina
```

## Example Workflows

### Scenario 1: First time run with all workflows

```bash
# Run everything
./run_all_workflows.sh --task-root /data/screening --proteins "PROTEIN1 PROTEIN2"

# Monitor progress
watch -n 60 ./run_all_workflows.sh status
```

### Scenario 2: Run only docking methods (skip structure prediction)

```bash
# Skip Boltz2, run only Vina and DiffDock
./run_all_workflows.sh --skip-boltz2 --skip-md-pbsa

# Or run them individually
./run_all_workflows.sh vina
./run_all_workflows.sh diffdock
```

### Scenario 3: Run only MD+PBSA (after DiffDock completes)

```bash
# Wait for DiffDock to complete, then run MD+PBSA
./run_all_workflows.sh md_pbsa
```

### Scenario 4: Multiple proteins with selective workflows

```bash
# Run Boltz2 and Vina for three proteins
./run_all_workflows.sh \
    --task-root /mnt/data \
    --proteins "JAK1 JAK2 JAK3" \
    --skip-diffdock \
    --skip-md-pbsa
```

## Troubleshooting

### "Script not found" error
Make sure all sub-scripts are in the same directory as the master script:
- `run_boltz2_batch.sh`
- `run_vina_batch.sh`
- `run_diffdock_batch.sh`
- `run_md_pbsa_batch.sh`

### "Script not executable" error
Make all scripts executable:
```bash
chmod +x run_all_workflows.sh
chmod +x run_boltz2_batch.sh
chmod +x run_vina_batch.sh
chmod +x run_diffdock_batch.sh
chmod +x run_md_pbsa_batch.sh
```

### Prerequisites check failed
Verify:
- Task root directory exists and is accessible
- Initial screening has been completed (`selected.csv` exists)
- Input files (PDB, GRO, TOP) exist for each protein

### No jobs submitted
Check:
- SLURM is available: `squeue`
- Partition names are correct in each sub-script
- Resource limits are within cluster quotas

## Integration with Existing Scripts

The master script works with your existing individual scripts. You can:

1. **Use master script for convenience**:
   ```bash
   ./run_all_workflows.sh
   ```

2. **Use individual scripts directly**:
   ```bash
   export MASTER_TASK_ROOT="/my/data"
   export MASTER_PROTEINS="PROTEIN1 PROTEIN2"
   ./run_boltz2_batch.sh
   ```

3. **Mix and match**:
   ```bash
   # Use master for some
   ./run_all_workflows.sh boltz2 vina

   # Use individual for others
   ./run_diffdock_batch.sh
   ```

All approaches work because sub-scripts check environment variables first, then fall back to defaults.
