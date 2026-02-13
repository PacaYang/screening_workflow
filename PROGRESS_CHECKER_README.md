# Screening Progress Checker

## Overview

The `check_screening_progress.sh` script provides a comprehensive overview of all screening workflows, showing completion status and failure rates for each algorithm.

## Features

- ✅ **AF3 (AlphaFold 3)** - Checks for `*_summary_confidences.json` files
- ✅ **Vina (AutoDock Vina)** - Checks for `docking_affinities.txt` files
- ✅ **DiffDock** - Checks for `rank1.sdf` files
- ✅ **MD+PBSA** - Checks both MD and PBSA phases separately
- ✅ **SLURM Jobs** - Shows running and pending jobs
- ✅ **Color-coded output** - Green (completed), Red (failed), Yellow (pending)
- ✅ **Success rates** - Percentage of successful completions
- ✅ **Sample outputs** - Shows example output files and values
- ✅ **Failed compound lists** - Lists up to 5 failed compounds per algorithm

## Usage

### Basic Usage

```bash
./check_screening_progress.sh
```

### With Custom Configuration

```bash
# Set task root and proteins via environment
export MASTER_TASK_ROOT="/path/to/data"
export MASTER_PROTEINS="PROTEIN1 PROTEIN2"
./check_screening_progress.sh
```

### Quick Check

```bash
# Quick one-liner to check overall progress
./check_screening_progress.sh | grep "Success rate"
```

## Output Format

### Summary Report

```
==========================================
  Screening Progress Report
==========================================
  Generated at: Mon Jan 15 10:00:00 UTC 2025
==========================================

==========================================
SUMMARY REPORT
==========================================
  Total proteins: 1
  Total compounds: 102

  Task root: /home/ubuntu/snake_test
  Proteins: JAK1JH1
```

### Per-Algorithm Status

```
==========================================
Protein: JAK1JH1
==========================================
  Total compounds in screening: 102

------------------------------------------
AlphaFold 3
------------------------------------------
  Total compounds: 102
  Completed: 95
  Failed: 3
  Pending: 4
  Success rate: 93%

  Sample output: JAK1JH1_0

  Failed compounds (invalid/empty JSON files):
    - JAK1JH1_12
    - JAK1JH1_45
    - JAK1JH1_78

------------------------------------------
AutoDock Vina
------------------------------------------
  Total compounds: 102 (across 10 parts)
  Completed: 98
  Failed: 2
  Pending: 2
  Success rate: 96%

  Sample output: input_part_0/lig5

------------------------------------------
DiffDock
------------------------------------------
  Total compounds: 102
  Completed: 90
  Failed: 5
  Pending: 7
  Success rate: 88%

  Sample output: JAK1JH1_0
  File size: 24K

  Failed compounds (missing or empty rank1.sdf):
    - JAK1JH1_8
    - JAK1JH1_23
    - JAK1JH1_47
    - JAK1JH1_68
    - JAK1JH1_89

------------------------------------------
MD + PBSA
------------------------------------------
  MD Phase:
    Completed: 85
    Failed: 2
    Pending: 15
    Success rate: 83%

  PBSA Phase:
    Completed: 80
    Failed: 3
    Pending: 19
    Success rate: 78%

  Overall (MD → PBSA pipeline):
    Total compounds: 102
    Fully completed (MD+PBSA): 80
    In progress (MD done, PBSA pending): 5

  Sample PBSA output: JAK1JH1_0
  Sample binding energy: -45.23 kcal/mol

------------------------------------------
SLURM Job Status
------------------------------------------
  AF3:       Running: 2, Pending: 2
  Vina:      Running: 0, Pending: 0
  DiffDock:  Running: 3, Pending: 4
  MD+PBSA:   Running: 8, Pending: 7

  Total:     Running: 13, Pending: 13
```

## Validation Criteria

### AF3
- ✅ **Success**: `*_summary_confidences.json` exists, is not empty, and contains "ranking_score"
- ❌ **Failed**: JSON file exists but is empty or missing "ranking_score"
- ⏳ **Pending**: JSON file doesn't exist

### Vina
- ✅ **Success**: `docking_affinities.txt` exists, is not empty, and contains "REMARK VINA RESULT"
- ❌ **Failed**: File exists but is empty or missing Vina results
- ⏳ **Pending**: File doesn't exist

### DiffDock
- ✅ **Success**: `rank1.sdf` exists and has content
- ❌ **Failed**: Directory exists but `rank1.sdf` is missing or empty
- ⏳ **Pending**: Directory doesn't exist

### MD Phase
- ✅ **Success**: `token.done` and `T298.xtc` exist with content
- ❌ **Failed**: Directory exists but missing token or trajectory
- ⏳ **Pending**: Directory doesn't exist

### PBSA Phase
- ✅ **Success**: `token.done` and `FINAL_RESULTS_MMPBSA.dat` exist with "DELTA TOTAL"
- ❌ **Failed**: Directory exists but missing token or results
- ⏳ **Pending**: Directory doesn't exist

## File Locations Checked

### AF3
```
${TASK_ROOT}/${PROTEIN}/fine_screening/AF3/output/${PROTEIN}_${i}/${PROTEIN}_${i}_summary_confidences.json
```

### Vina
```
${TASK_ROOT}/${PROTEIN}/fine_screening/Vina/output/input_part_${i}/lig${j}/docking_affinities.txt
```

### DiffDock
```
${TASK_ROOT}/${PROTEIN}/fine_screening/PBSA/DiffDock/output/${PROTEIN}_${i}/rank1.sdf
```

### MD
```
${TASK_ROOT}/${PROTEIN}/fine_screening/PBSA/PBSA/MD/${PROTEIN}_${i}/token.done
${TASK_ROOT}/${PROTEIN}/fine_screening/PBSA/PBSA/MD/${PROTEIN}_${i}/T298.xtc
```

### PBSA
```
${TASK_ROOT}/${PROTEIN}/fine_screening/PBSA/PBSA/PBSA/${PROTEIN}_${i}/token.done
${TASK_ROOT}/${PROTEIN}/fine_screening/PBSA/PBSA/PBSA/${PROTEIN}_${i}/FINAL_RESULTS_MMPBSA.dat
```

## Use Cases

### 1. Daily Progress Check
```bash
# Run once per day to monitor overall progress
./check_screening_progress.sh > daily_report_$(date +%Y%m%d).txt
```

### 2. Find Failed Compounds
```bash
# Extract list of failed compounds
./check_screening_progress.sh | grep -A 10 "Failed compounds"
```

### 3. Check Success Rates
```bash
# See only success rates
./check_screening_progress.sh | grep "Success rate"
```

### 4. Monitor Active Jobs
```bash
# See only SLURM job status
./check_screening_progress.sh | grep -A 10 "SLURM Job Status"
```

### 5. Compare Multiple Proteins
```bash
# Check progress for multiple proteins
export MASTER_PROTEINS="JAK1 JAK2 JAK3"
./check_screening_progress.sh
```

### 6. Continuous Monitoring
```bash
# Watch progress every 5 minutes
watch -n 300 ./check_screening_progress.sh
```

## Integration with Other Scripts

This script complements the individual monitoring scripts:

```bash
# Overall progress across all algorithms
./check_screening_progress.sh

# Detailed monitoring for specific workflows
./monitor_af3_jobs.sh      # AF3 details
./monitor_vina_jobs.sh     # Vina details
./monitor_diffdock_jobs.sh # DiffDock details
./monitor_md_pbsa_jobs.sh  # MD+PBSA details

# Master control and status
./run_all_workflows.sh status
```

## Interpreting Results

### High Success Rate (>90%)
✅ Workflow is running smoothly
- Continue monitoring
- No action needed

### Medium Success Rate (70-90%)
⚠️ Some failures occurring
- Check error logs for failed compounds
- Investigate common failure patterns
- May need to adjust parameters

### Low Success Rate (<70%)
❌ Significant issues
- Check SLURM error logs
- Verify input files
- Check resource availability (disk space, memory)
- Review workflow configuration

### High Pending Count
⏳ Jobs are queued or running
- Check SLURM queue: `squeue -u $USER`
- Monitor with: `watch -n 60 'squeue -u $USER'`
- Wait for jobs to complete

## Troubleshooting

### No output directories found
```
Status: Not started (output directory not found)
```
**Solution**: Workflow hasn't been submitted yet. Run the appropriate submission script.

### All compounds showing as pending
```
Completed: 0
Failed: 0
Pending: 102
```
**Solution**: Jobs are still running or queued. Check SLURM status section.

### High failure rate
```
Completed: 20
Failed: 80
Pending: 2
Success rate: 20%
```
**Solution**:
1. Check error logs in output directories
2. Look for common errors (out of memory, disk space, etc.)
3. Verify input file quality
4. Check if jobs are timing out

### Script shows error about selected.csv
```
ERROR: Selected compounds file not found
```
**Solution**: Initial screening hasn't been completed. Run initial screening first.

## Color Legend

When viewing in a color-enabled terminal:
- 🟢 **Green** - Completed successfully
- 🔴 **Red** - Failed
- 🟡 **Yellow** - Pending/In progress

## Performance

The script is optimized for fast execution:
- Uses `find` with limits to avoid scanning too many files
- Checks only necessary files
- Shows sample outputs instead of processing all
- Typical runtime: <30 seconds for 100 compounds per protein

## Example Output Interpretation

```
MD + PBSA
  MD Phase:
    Completed: 85
    Pending: 15

  PBSA Phase:
    Completed: 80
    Pending: 19
```

**Interpretation:**
- 85 compounds completed MD
- 80 compounds completed both MD and PBSA
- 5 compounds finished MD but waiting for PBSA (85 - 80)
- 15 compounds haven't started MD yet
- 19 compounds haven't completed PBSA (includes the 15 not started + 4 others)
