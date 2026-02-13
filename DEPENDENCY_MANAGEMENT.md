# Workflow Dependency Management

## Overview

The master script (`run_all_workflows.sh`) now supports automatic dependency management between DiffDock and MD+PBSA workflows.

## Feature: Wait for DiffDock Before MD+PBSA

### Configuration

Edit lines 35-36 in `run_all_workflows.sh`:

```bash
# Dependency control - wait for DiffDock before MD+PBSA
WAIT_FOR_DIFFDOCK=1  # Set to 1 to wait for DiffDock completion before running MD+PBSA
WAIT_POLL_INTERVAL=300  # Check every 5 minutes (300 seconds)
```

### How It Works

When `WAIT_FOR_DIFFDOCK=1`:

1. **If running `all` command with DiffDock enabled:**
   - Submits Boltz2 jobs (if enabled)
   - Submits Vina jobs (if enabled)
   - Submits DiffDock jobs
   - **Waits for all DiffDock jobs to complete**
   - Submits MD+PBSA jobs after DiffDock finishes

2. **If running `all` command with DiffDock disabled:**
   - Checks if DiffDock is already complete
   - If complete: runs MD+PBSA
   - If not complete: skips MD+PBSA with warning

3. **If running individual workflow commands** (e.g., `./run_all_workflows.sh md_pbsa`):
   - Ignores the wait setting
   - Runs immediately

When `WAIT_FOR_DIFFDOCK=0`:
- MD+PBSA runs immediately without checking DiffDock status
- Use this if you want to submit all jobs at once (MD+PBSA jobs will skip if SDF files are missing)

## Usage Examples

### Example 1: Automatic Dependency Management (Default)

```bash
# Run all workflows - MD+PBSA waits for DiffDock
./run_all_workflows.sh

# Output:
# [Submits Boltz2 jobs]
# [Submits Vina jobs]
# [Submits DiffDock jobs]
#
# ==========================================
# MD+PBSA has dependency on DiffDock
# ==========================================
# Checking DiffDock completion status...
# JAK1JH1: DiffDock 5/10 parts completed
#
# DiffDock is not complete. Waiting for completion...
# Check interval: 300 seconds
# Press Ctrl+C to stop waiting and skip MD+PBSA
#
# [Waits 5 minutes]
#
# ==========================================
# Wait check #1 at Mon Jan 15 10:05:00 UTC 2024
# ==========================================
# JAK1JH1: DiffDock 8/10 parts completed
#
# Still waiting... (checked 1 times)
# Next check in 300 seconds
#
# [Continues checking until complete]
#
# ==========================================
# DiffDock completed! Proceeding to MD+PBSA
# ==========================================
# [Submits MD+PBSA jobs]
```

### Example 2: Skip Waiting (Submit All at Once)

Edit script to set `WAIT_FOR_DIFFDOCK=0`:

```bash
# Run all workflows - all jobs submitted immediately
./run_all_workflows.sh

# MD+PBSA jobs will be submitted immediately
# Jobs with missing SDF files will skip gracefully when they run
```

### Example 3: Run MD+PBSA After DiffDock (Without Running DiffDock)

```bash
# DiffDock already completed in a previous run
./run_all_workflows.sh --skip-boltz2 --skip-vina --skip-diffdock

# Output:
# Skipping Boltz2 (disabled)
# Skipping Vina (disabled)
# Skipping DiffDock (disabled)
#
# ==========================================
# Checking DiffDock completion before MD+PBSA
# ==========================================
# JAK1JH1: DiffDock 10/10 parts completed
#
# ==========================================
# DiffDock is already complete for all proteins
# ==========================================
# [Submits MD+PBSA jobs]
```

### Example 4: Manual Control with Individual Commands

```bash
# Submit DiffDock
./run_all_workflows.sh diffdock

# Wait manually (monitor progress)
./monitor_diffdock_jobs.sh

# When complete, submit MD+PBSA
./run_all_workflows.sh md_pbsa
```

### Example 5: Interrupt Waiting

```bash
# Start waiting for DiffDock
./run_all_workflows.sh

# While waiting, press Ctrl+C to stop

# Output:
# Still waiting... (checked 3 times)
# Next check in 300 seconds
# ^C
#
# ==========================================
# DiffDock wait was interrupted. MD+PBSA not started.
# ==========================================
```

## Monitoring During Wait

While the script is waiting for DiffDock, you can monitor progress in another terminal:

```bash
# Terminal 1: Master script waiting
./run_all_workflows.sh

# Terminal 2: Monitor DiffDock progress
./monitor_diffdock_jobs.sh

# Or check SLURM queue
squeue -u $USER | grep diffdock
```

## Completion Detection

The script checks if DiffDock is complete by:

1. Counting expected parts from `${TASK_ROOT}/${PROTEIN}/fine_screening/Vina/input/*.csv`
2. Counting completed parts from `${TASK_ROOT}/${PROTEIN}/fine_screening/PBSA/DiffDock/output/*.done` files
3. Comparing: If `completed >= expected` for all proteins, DiffDock is considered complete

## Benefits

### With `WAIT_FOR_DIFFDOCK=1` (Recommended)
✅ Ensures MD+PBSA only runs when DiffDock is complete
✅ Maximizes success rate (no missing SDF files)
✅ Single command to run entire pipeline
✅ Automatic dependency management
✅ Can interrupt and resume later

### With `WAIT_FOR_DIFFDOCK=0`
✅ Submit all jobs at once
✅ Faster submission (no waiting)
⚠️ MD+PBSA jobs may skip if DiffDock not done
⚠️ Need to resubmit MD+PBSA for failed compounds

## Adjusting Wait Interval

Modify `WAIT_POLL_INTERVAL` based on expected DiffDock runtime:

```bash
# Check every 2 minutes (for fast testing)
WAIT_POLL_INTERVAL=120

# Check every 5 minutes (default, balanced)
WAIT_POLL_INTERVAL=300

# Check every 15 minutes (for long-running jobs)
WAIT_POLL_INTERVAL=900

# Check every hour (for very long workflows)
WAIT_POLL_INTERVAL=3600
```

## Troubleshooting

### Script never detects completion
- Check that DiffDock jobs are creating `.done` token files
- Verify paths are correct
- Run: `ls ${TASK_ROOT}/${PROTEIN}/fine_screening/PBSA/DiffDock/output/*.done`

### False completion detected
- Verify input CSV files exist
- Check if Vina/input directory has the expected CSV files

### Want to skip waiting and force MD+PBSA
1. Press Ctrl+C to stop waiting
2. Set `WAIT_FOR_DIFFDOCK=0` in the script
3. Run: `./run_all_workflows.sh md_pbsa`

## Implementation Details

The dependency check is implemented through three functions:

1. **`check_diffdock_complete()`** - Checks if all DiffDock jobs are done
2. **`wait_for_diffdock()`** - Polls until DiffDock completes
3. **Main workflow logic** - Calls wait function before MD+PBSA

The wait function uses `sleep` and a loop to periodically check completion status.
