#!/bin/bash
# Example usage of the peptide design pipeline

# Example 1: Basic usage with hotspot residues
bash run_peptide_pipeline.sh \
  --protein JAK1 \
  --protein-pdb /path/to/jak1.pdb \
  --peptide-length 15 \
  --output-dir /home/yangl_pacagen_com/peptide_results/jak1_peptides \
  --hotspot-residues "A:45,A:67,A:123"

# Example 2: Quick test run (fewer designs for testing)
bash run_peptide_pipeline.sh \
  --protein JAK1 \
  --protein-pdb /path/to/jak1.pdb \
  --peptide-length 10 \
  --output-dir /home/yangl_pacagen_com/test_peptides \
  --n-designs 5 \
  --n-seqs 2 \
  --n-batches 2

# Example 3: No hotspot constraints (explore full surface)
bash run_peptide_pipeline.sh \
  --protein TARGET \
  --protein-pdb /path/to/target.pdb \
  --peptide-length 12 \
  --output-dir /path/to/output

# Example 4: Resume from streaming collection (step 5)
bash run_peptide_pipeline.sh \
  --protein JAK1 \
  --protein-pdb /path/to/jak1.pdb \
  --peptide-length 15 \
  --output-dir /home/yangl_pacagen_com/peptide_results/jak1_peptides \
  --resume-from 5

# Monitoring commands
# -------------------

# Check running SLURM jobs
squeue -u $USER | grep af3

# Watch streaming results
watch -n 60 'ls -lh /path/to/output/05_results/streaming_updates/'

# Check latest batch results
tail -20 /path/to/output/05_results/streaming_updates/batch_*.csv | sort -t, -k8 -rn

# View final rankings
head -20 /path/to/output/05_results/ranked_peptides.csv

# Check pipeline state
cat /path/to/output/.pipeline_state.json | jq .

# Troubleshooting
# ---------------

# View RFDiffusion logs
tail -f /path/to/output/01_rfdiffusion/logs/rfdiffusion.log

# View ProteinMPNN logs
tail -f /path/to/output/02_proteinmpnn/logs/proteinmpnn.log

# Check AF3 job logs
ls /path/to/output/04_af3_output/logs/

# View collection state
cat /path/to/output/04_af3_output/.collection_state/AF3_peptide.json | jq .
