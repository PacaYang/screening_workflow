# Nesso, Protenix, and OpenFold3 score collection

The three collectors read native prediction JSONs and write one `summary.csv`
row per row of `initial_screening/selected.csv`. They use the selected row index
as the prediction key, so repeated SMILES remain separate predictions.

```bash
python scripts/scoring/nesso_scores.py \
  --results-folder /path/to/context/fine_screening/Nesso/output \
  --selected-csv /path/to/context/initial_screening/selected.csv \
  --output-dir /path/to/context/fine_screening/Nesso

python scripts/scoring/protenix_scores.py \
  --results-folder /path/to/context/fine_screening/Protenix/output \
  --selected-csv /path/to/context/initial_screening/selected.csv \
  --output-dir /path/to/context/fine_screening/Protenix

python scripts/scoring/openfold3_scores.py \
  --results-folder /path/to/context/fine_screening/OpenFold3/output \
  --selected-csv /path/to/context/initial_screening/selected.csv \
  --output-dir /path/to/context/fine_screening/OpenFold3
```

Protenix selects the highest `ranking_score` across valid samples. OpenFold3
selects the highest `sample_ranking_score`. The chosen seed, sample number,
valid sample count, and relative source JSON path are recorded. Ties select
the lowest seed and sample number. Nesso copies its reported affinity,
binding probability, and entropy values; no score conversion is applied.

Rows with no valid native score have blank score columns and a `status` and
`reason`. `missing_score` means a job was marked successful in the ledger but
produced no score JSON. `failed_job` comes from the job ledger. `pending`
means no terminal job result is recorded. `unsupported_group` means the input
manifest assigned multiple ligands to one model job; this collector does not
reuse a complex score for individual ligands. Collectors rebuild the CSV on
each invocation so retries update existing rows.

To collect a benchmark run with `predictions.csv` and multiple `ctx_*`
directories, run:

```bash
python scripts/scoring/collect_ml_run_scores.py --run-root /path/to/run
```

This writes three method summaries under every context and one
`/path/to/run/compiled_scores.csv` with prediction metadata and prefixed
scores. The run-wide join uses context and compound index. For a standard
pipeline target, Stage 6 calls the collectors and `compile_scores.py` writes
`fine_screening/compiled_scores.csv` from whichever enabled methods have
summaries. `compile_scores.py` joins by selected row index when
`selected.csv` exists; without it, SMILES must be unique in every input.
