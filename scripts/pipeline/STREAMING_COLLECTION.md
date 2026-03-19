# Streaming Results Collection

## Overview

Stage 6 now supports **streaming results collection**, which monitors and collects results incrementally as fine screening jobs complete, rather than waiting for all jobs to finish.

## Benefits

- **Faster feedback** - See results as they complete, not after all jobs finish
- **Better resource utilization** - Results collected while jobs still running
- **Fault tolerance** - Resume collection after failures without reprocessing
- **Early problem detection** - Warnings about empty/incomplete results during run
- **Progress visibility** - Real-time tracking of collection progress
- **Reduced wait time** - No need to wait for slowest job before seeing any results

## Usage

### Streaming Mode (Default)

```bash
# Use streaming collection with default 1-hour interval
bash scripts/pipeline/run_full_pipeline.sh \
  --task-root /path/to/task \
  --proteins "PROTEIN1 PROTEIN2" \
  --start-from 6 \
  run

# Custom collection interval (30 minutes)
bash scripts/pipeline/run_full_pipeline.sh \
  --task-root /path/to/task \
  --proteins "PROTEIN1" \
  --collection-interval 1800 \
  run

# Limit to 5 collection iterations
bash scripts/pipeline/run_full_pipeline.sh \
  --task-root /path/to/task \
  --proteins "PROTEIN1" \
  --max-collection-iter 5 \
  run
```

### Batch Mode (Original Behavior)

```bash
# Use old batch collection (wait for all jobs)
bash scripts/pipeline/run_full_pipeline.sh \
  --task-root /path/to/task \
  --proteins "PROTEIN1" \
  --batch-collection \
  run
```

## Configuration Options

| Option | Default | Description |
|--------|---------|-------------|
| `--streaming-collection` | enabled | Use streaming collection mode |
| `--batch-collection` | disabled | Use batch collection mode (original) |
| `--collection-interval N` | 3600 | Seconds between collection runs |
| `--max-collection-iter N` | 0 | Max iterations (0=infinite) |

## How It Works

### Collection State Tracking

For each method (AF3, Boltz2, Vina, PBSA), the system maintains a state file:

```
${TASK_ROOT}/${PROTEIN}/fine_screening/.collection_state/${METHOD}.json
```

Example state file:
```json
{
  "method": "AF3",
  "protein": "JAK1JH1",
  "last_collection": "2026-03-06T22:30:00Z",
  "total_expected": 3000,
  "collected_compounds": {
    "0": {
      "timestamp": "2026-03-06T20:15:00Z",
      "status": "collected",
      "has_data": true
    },
    "1": {
      "timestamp": "2026-03-06T20:16:00Z",
      "status": "collected",
      "has_data": false,
      "warning": "Missing required JSON files"
    }
  },
  "summary_csv_rows": 1
}
```

### Collection Loop

1. **Check for new results** - Scan output directories for completed compounds
2. **Filter uncollected** - Compare against state file to find new results
3. **Extract scores** - Call scoring scripts in incremental mode
4. **Append to CSV** - Add new rows to existing summary files
5. **Update state** - Record which compounds have been collected
6. **Validate** - Check for empty or incomplete results
7. **Wait** - Sleep for configured interval before next iteration
8. **Repeat** - Continue until all compounds collected or max iterations reached

### Incremental Scoring Scripts

The scoring scripts (`af3_scores.py`, `boltz2_scores.py`) now support:

- `--incremental` - Enable incremental collection mode
- `--state-file` - Path to collection state JSON
- `--append` - Append to existing CSV instead of overwriting

Example:
```bash
python scripts/scoring/af3_scores.py \
  --af3-results-folder /path/to/output \
  --output-dir /path/to/results \
  --incremental \
  --state-file /path/to/state.json \
  --append
```

## Result Validation

The system automatically validates collected results:

- **Empty CSV detection** - Warns if summary CSV has 0 rows
- **Low collection rate** - Warns if less than 50% of expected compounds collected
- **Missing data** - Tracks compounds with missing or invalid files
- **Progress reporting** - Shows X/Y compounds collected

Example output:
```
[INFO] JAK1JH1/AF3: Collecting 150 uncollected compounds (2850/3000 done)
[INFO] JAK1JH1/AF3: 2900/3000 compounds collected
[WARN] JAK1JH1/Boltz2: Only 1200/3000 compounds collected (40%)
```

## File Structure

```
scripts/
├── pipeline/
│   ├── lib/
│   │   └── collection_state.sh          # State management functions
│   ├── stages/
│   │   ├── stage6_collect_results.sh    # Main entry point (delegates to streaming/batch)
│   │   └── stage6_collect_results_streaming.sh  # Streaming collection implementation
│   └── run_full_pipeline.sh             # Master controller
└── scoring/
    ├── af3_scores.py                    # AF3 scoring (incremental support)
    └── boltz2_scores.py                 # Boltz2 scoring (incremental support)
```

## Implementation Status

| Method | Incremental Support | Status |
|--------|-------------------|--------|
| AF3 | ✅ Yes | Implemented |
| Boltz2 | ✅ Yes | Implemented |
| Vina | ⏳ Pending | Not yet implemented |
| PBSA | ⏳ Pending | Not yet implemented |

## Troubleshooting

### Collection not progressing

Check if jobs are actually completing:
```bash
# Check job status
squeue -u $USER

# Check output directories
ls ${TASK_ROOT}/${PROTEIN}/fine_screening/AF3/output/
```

### State file corruption

Reset state and re-collect:
```bash
# Remove state files
rm -rf ${TASK_ROOT}/${PROTEIN}/fine_screening/.collection_state/

# Remove existing summary (will be regenerated)
rm ${TASK_ROOT}/${PROTEIN}/fine_screening/AF3/summary.csv

# Re-run collection
bash scripts/pipeline/stages/stage6_collect_results_streaming.sh \
  --task-root ${TASK_ROOT} \
  --proteins ${PROTEIN} \
  --max-iterations 1
```

### Duplicate rows in CSV

This shouldn't happen with proper state tracking, but if it does:
```bash
# Remove duplicates
python -c "
import pandas as pd
df = pd.read_csv('summary.csv')
df = df.drop_duplicates(subset=['SMILES', 'folder'])
df.to_csv('summary_dedup.csv', index=False)
"
```

## Performance Considerations

- **Collection interval**: Shorter intervals = more frequent updates but more overhead
- **Max iterations**: Set to limit total runtime (e.g., 24 iterations = 24 hours with 1-hour interval)
- **Concurrent safety**: Multiple collection processes can run safely (state file uses atomic updates)

## Migration from Batch Mode

Existing pipelines will automatically use streaming mode. To opt out:

```bash
# Add to your pipeline invocation
--batch-collection
```

Or set in configuration file:
```yaml
collection:
  mode: batch  # or "streaming"
  interval: 3600
```
