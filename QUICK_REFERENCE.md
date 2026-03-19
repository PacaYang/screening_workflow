# Quick Reference: Streaming Results Collection

## Quick Start

### Use streaming collection (default)
```bash
bash scripts/pipeline/run_full_pipeline.sh \
  --task-root /path/to/task \
  --proteins "PROTEIN1" \
  run
```

### Use batch collection (old behavior)
```bash
bash scripts/pipeline/run_full_pipeline.sh \
  --task-root /path/to/task \
  --proteins "PROTEIN1" \
  --batch-collection \
  run
```

## Common Options

| Option | Default | Description |
|--------|---------|-------------|
| `--streaming-collection` | enabled | Use streaming mode |
| `--batch-collection` | disabled | Use batch mode |
| `--collection-interval 1800` | 3600 | Check every 30 min |
| `--max-collection-iter 10` | 0 | Stop after 10 iterations |

## Check Progress

```bash
# View collection state
cat ${TASK_ROOT}/${PROTEIN}/fine_screening/.collection_state/AF3.json | jq .

# Count collected compounds
jq '.collected_compounds | length' ${TASK_ROOT}/${PROTEIN}/fine_screening/.collection_state/AF3.json

# Check for warnings
jq '.collected_compounds | to_entries[] | select(.value.has_data == false)' \
  ${TASK_ROOT}/${PROTEIN}/fine_screening/.collection_state/AF3.json

# View summary CSV
head ${TASK_ROOT}/${PROTEIN}/fine_screening/AF3/summary.csv
wc -l ${TASK_ROOT}/${PROTEIN}/fine_screening/AF3/summary.csv
```

## Troubleshooting

### Reset collection state
```bash
# Remove state files
rm -rf ${TASK_ROOT}/${PROTEIN}/fine_screening/.collection_state/

# Remove summary CSV (will be regenerated)
rm ${TASK_ROOT}/${PROTEIN}/fine_screening/AF3/summary.csv

# Re-run collection
bash scripts/pipeline/stages/stage6_collect_results_streaming.sh \
  --task-root ${TASK_ROOT} \
  --proteins ${PROTEIN} \
  --max-iterations 1
```

### Force batch collection once
```bash
bash scripts/pipeline/stages/stage6_collect_results.sh \
  --task-root ${TASK_ROOT} \
  --proteins ${PROTEIN} \
  --batch-collection
```

## Files and Locations

```
${TASK_ROOT}/${PROTEIN}/fine_screening/
├── .collection_state/          # State tracking
│   ├── AF3.json               # AF3 collection state
│   ├── Boltz2.json            # Boltz2 collection state
│   ├── Vina.json              # Vina collection state (future)
│   └── PBSA.json              # PBSA collection state (future)
├── AF3/
│   ├── output/                # Raw results
│   │   ├── 0/                 # Compound folders
│   │   ├── 1/
│   │   └── ...
│   └── summary.csv            # Collected results (appended incrementally)
├── Boltz2/
│   ├── output/
│   └── summary.csv
├── Vina/
│   ├── output/
│   └── results.csv
└── PBSA/
    ├── PBSA/PBSA/
    └── summary.csv
```

## Implementation Status

| Method | Incremental Support |
|--------|-------------------|
| AF3 | ✅ Implemented |
| Boltz2 | ✅ Implemented |
| Vina | ⏳ Pending |
| PBSA | ⏳ Pending |

## Key Scripts

- `scripts/pipeline/lib/collection_state.sh` - State management
- `scripts/pipeline/stages/stage6_collect_results_streaming.sh` - Streaming collection
- `scripts/pipeline/stages/stage6_collect_results.sh` - Entry point (delegates)
- `scripts/scoring/af3_scores.py` - AF3 incremental scoring
- `scripts/scoring/boltz2_scores.py` - Boltz2 incremental scoring

## Documentation

- Full guide: `scripts/pipeline/STREAMING_COLLECTION.md`
- Implementation summary: `IMPLEMENTATION_SUMMARY.md`
- Test script: `test_streaming_collection.sh`
