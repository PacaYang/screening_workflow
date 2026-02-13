# Boltz2 Batch Prediction Scripts

This directory contains automation scripts for running Boltz2 predictions in batch mode using SLURM job scheduling.

## Overview

The workflow splits Boltz2 predictions into **40 batches** by default, allowing parallel processing on HPC clusters. Each batch processes multiple SMILES compounds from the initial screening results.

## Files

- `run_boltz2_batch.sh` - Main submission script that creates and submits SLURM jobs
- `monitor_boltz2_jobs.sh` - Monitoring script to track job progress

## Prerequisites

Before running the prediction step, ensure the following are completed:

1. **Initial screening** - The initial screening must be completed with results in:
   - `${TASK_ROOT}/${PROTEIN}/initial_screening/selected.csv`

2. **Input preparation** - Boltz2 input YAML files must be generated:
   - Token file: `/home/ubuntu/${PROTEIN}/boltz2_tmp/boltz_input.done`
   - Input directory: `/home/ubuntu/${PROTEIN}/boltz2_tmp/input/`

3. **Prefold completed** - Protein prefolding with MSA must be done first

## Configuration

Edit `run_boltz2_batch.sh` to customize:

```bash
# SLURM settings
PARTITION="gpu"           # Your GPU partition name
TIME_LIMIT="04:00:00"     # Time limit per batch job
MEMORY="32G"              # Memory per job
CPUS_PER_TASK=4          # Number of CPUs
GPUS_PER_NODE=1          # Number of GPUs per job

# Batch settings
N_BATCHES=40             # Number of batches to split work into
```

Update the protein list in the `get_proteins()` function:
```bash
get_proteins() {
    echo "JAK1JH1 PROTEIN2 PROTEIN3"  # Space-separated list
}
```

## Usage

### 1. Submit Batch Jobs

```bash
./run_boltz2_batch.sh
```

This will:
- Calculate the number of SMILES compounds to process
- Split work into 40 batches
- Submit one SLURM job per batch
- Save job IDs to `${TASK_ROOT}/${PROTEIN}/fine_screening/Boltz2/output/submitted_jobs.txt`

### 2. Monitor Progress

```bash
./monitor_boltz2_jobs.sh
```

This displays:
- Number of batches submitted vs completed
- SLURM job status (running, pending)
- Recent completions
- Progress percentage
- Warning if errors detected in logs

You can also use standard SLURM commands:
```bash
# View all your jobs
squeue -u $USER

# View specific protein's jobs
squeue -u $USER | grep boltz2_JAK1JH1

# Cancel a specific job
scancel <JOB_ID>

# Cancel all Boltz2 jobs for a protein
scancel -n boltz2_JAK1JH1_*
```

## Output Structure

```
${TASK_ROOT}/${PROTEIN}/fine_screening/Boltz2/output/
├── batch_0.tar.gz              # Compressed results for batch 0
├── batch_1.tar.gz              # Compressed results for batch 1
├── ...
├── batch_39.tar.gz             # Compressed results for batch 39
├── slurm_batch_0_<jobid>.out  # SLURM stdout for batch 0
├── slurm_batch_0_<jobid>.err  # SLURM stderr for batch 0
├── ...
├── submitted_jobs.txt          # List of submitted job IDs
└── token/
    ├── batch_0.done           # Completion markers
    ├── batch_1.done
    └── ...
```

Each tar.gz file contains:
- `affinity_*.json` - Affinity predictions
- `confidence*.json` - Confidence scores
- `*model*.cif` - Structure files

## Workflow Integration

This script corresponds to the **boltz2** rule in the Snakefile (lines 643-718). After all batches complete:

1. The Snakemake workflow can continue with the `collect_boltz2` rule
2. Results will be decompressed and analyzed to create `summary.csv`

## Resource Optimization

The script uses local temporary storage (`/tmp`) on compute nodes to:
- Reduce I/O load on shared filesystems
- Speed up intermediate file operations
- Only transfer final compressed results to shared storage

This is especially important for Boltz2, which generates many intermediate files.

## Troubleshooting

### Jobs not starting?
- Check GPU availability: `sinfo -p gpu`
- Verify partition name in the script matches your cluster

### Jobs failing?
- Check error logs: `${OUTPUT_DIR}/slurm_batch_*_*.err`
- Verify conda environment: `conda activate boltz`
- Ensure input YAML files exist: `ls /home/ubuntu/${PROTEIN}/boltz2_tmp/input/`

### Out of memory?
- Increase `MEMORY` parameter in the script
- Consider reducing batch size (increase `N_BATCHES`)

### Missing results?
- Check if jobs completed successfully
- Verify token files exist: `ls ${TOKEN_DIR}/`
- Look for tar.gz files: `ls ${OUTPUT_DIR}/batch_*.tar.gz`

## Example Run

```bash
# 1. Submit jobs
$ ./run_boltz2_batch.sh
[2025-11-24 10:00:00] INFO: Starting Boltz2 batch job submission
[2025-11-24 10:00:00] INFO: Number of batches per protein: 40
[2025-11-24 10:00:00] INFO: Processing protein: JAK1JH1
[2025-11-24 10:00:00] INFO: Total SMILES for JAK1JH1: 102
[2025-11-24 10:00:00] INFO: Batch size: 3
[2025-11-24 10:00:01] INFO: Submitted batch 0 for protein JAK1JH1 (Job ID: 12345)
[2025-11-24 10:00:02] INFO: Submitted batch 1 for protein JAK1JH1 (Job ID: 12346)
...

# 2. Monitor progress
$ ./monitor_boltz2_jobs.sh
========================================
Boltz2 Batch Job Monitoring
========================================

Protein: JAK1JH1
----------------------------------------
Total batches submitted: 40
Completed batches: 15
Compressed outputs: 15
Progress: 37%

SLURM Job Status:
  Running: 8
  Pending: 17

Recently completed batches (last 5):
  - batch_14.done
  - batch_13.done
  - batch_12.done
  - batch_11.done
  - batch_10.done
```

## Notes

- Each batch processes approximately `N_SMILES / 40` compounds
- Adjust `N_BATCHES` based on your cluster's GPU availability and job limits
- The script uses `--override` flag to rerun predictions if needed
- Compressed archives save ~90% storage space vs uncompressed outputs
