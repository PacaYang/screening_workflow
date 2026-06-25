# Peptide Design Pipeline

RFDiffusion → ProteinMPNN → AF3 validation pipeline for designing peptide binders.

## Quick Start

```bash
bash run_peptide_pipeline.sh \
  --protein JAK1 \
  --protein-pdb /path/to/jak1.pdb \
  --peptide-length 15 \
  --output-dir /path/to/output \
  --hotspot-residues "A:45,A:67,A:123"
```

## Pipeline Overview

1. **RFDiffusion** (`step1_rfdiffusion.sh`): Generates peptide backbones that bind to the target protein
2. **ProteinMPNN** (`step2_proteinmpnn.sh`): Designs sequences for generated backbones
3. **Prepare AF3** (`step3_prepare_af3.py`): Creates AF3 JSON input files
4. **Submit AF3** (`step4_submit_af3.sh`): Submits SLURM jobs for structure prediction
5. **Stream Collection** (`step5_collect_streaming.sh`): Monitors and collects results incrementally
6. **Finalize** (`step6_finalize_ranking.py`): Ranks peptides by AF3 confidence scores

## Parameters

**Required:**
- `--protein`: Protein name
- `--protein-pdb`: Path to protein PDB structure
- `--peptide-length`: Target peptide length (e.g., 10-20)
- `--output-dir`: Output directory

**Optional:**
- `--hotspot-residues`: Binding site residues (e.g., "A:45,A:67")
- `--n-designs`: Number of backbones (default: 50)
- `--n-seqs`: Sequences per backbone (default: 10)
- `--n-batches`: AF3 batch count (default: 40)
- `--resume-from`: Resume from step number (1-6)

## Output Structure

```
${OUTPUT_DIR}/
├── 01_rfdiffusion/backbones/    # Generated PDB backbones
├── 02_proteinmpnn/sequences/    # Designed sequences
├── 03_af3_input/                # AF3 JSON inputs
├── 04_af3_output/predictions/   # AF3 predictions
└── 05_results/
    ├── ranked_peptides.csv      # Final ranked results
    └── streaming_updates/       # Incremental results
```

## Ranking Criteria

Peptides are ranked by: **0.4 × ipTM + 0.3 × pTM + 0.3 × pLDDT**

- **ipTM**: Interface confidence (binding quality)
- **pTM**: Overall structure confidence
- **pLDDT**: Per-residue confidence

## Resume Capability

Pipeline tracks progress in `.pipeline_state.json`. Resume after interruption:

```bash
bash run_peptide_pipeline.sh --resume-from 5 [other args...]
```

## Monitoring

- Check SLURM jobs: `squeue -u $USER | grep af3`
- Watch streaming results: `tail -f ${OUTPUT_DIR}/05_results/streaming_updates/batch_*.csv`
- Check logs: `${OUTPUT_DIR}/04_af3_output/logs/`

## Requirements

- RFDiffusion installed in `~/Applications/RFdiffusion`
- ProteinMPNN installed in `~/Applications/ProteinMPNN`
- AlphaFold3 installed in `~/Applications/alphafold3`
- Conda environments: `SE3nv` (RFDiffusion), `mlfold` (ProteinMPNN), `af3_test` (AF3)
