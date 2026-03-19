# Binding Site Analysis for Screening Workflow

This implementation adds binding site analysis capabilities to the existing scoring scripts for AF3, Boltz2, and RoseTTAFold predictions.

## Overview

The implementation consists of two phases:

### Phase 1: Enhanced Scoring Scripts
Modified existing scoring scripts to extract binding site information:
- `af3_scores.py` - Extracts binding sites from AF3 CIF files
- `boltz2_scores.py` - Extracts binding sites from Boltz2 CIF files
- `rosettafold_scores.py` - Extracts binding sites from RoseTTAFold PDB files

Each script now adds three new columns to the summary.csv output:
- `binding_site_residues`: Comma-separated list of residues within 5Å of ligand (e.g., "A:CYS45,A:ARG67")
- `binding_site_center`: X,Y,Z coordinates of ligand center of mass
- `num_binding_residues`: Count of residues in binding site

### Phase 2: Binding Site Analysis Scripts
Created three new analysis scripts that perform clustering and visualization:
- `analyze_af3_binding_sites.py`
- `analyze_boltz2_binding_sites.py`
- `analyze_rosettafold_binding_sites.py`

## Installation

Install required dependencies:

```bash
pip install -r scripts/scoring/binding_site_requirements.txt
```

Required packages:
- biopython (for structure parsing)
- scikit-learn (for DBSCAN clustering)
- scipy (for distance calculations)
- matplotlib (for visualizations)
- seaborn (for heatmaps)

## Usage

### Step 1: Generate Enhanced Summary Files

Run the existing scoring scripts as usual. They will now automatically extract binding site information:

```bash
# AF3
python scripts/scoring/af3_scores.py \
  --af3-results-folder ~/Projects/B3/IL6/IL6RB/fine_screening/AF3/results \
  --output-dir ~/Projects/B3/IL6/IL6RB/fine_screening/AF3

# Boltz2
python scripts/scoring/boltz2_scores.py \
  --boltz-results-folder ~/Projects/B3/IL6/IL6RB/fine_screening/Boltz2/output \
  --output-dir ~/Projects/B3/IL6/IL6RB/fine_screening/Boltz2

# RoseTTAFold
python scripts/scoring/rosettafold_scores.py \
  --rfaa-results-folder ~/Projects/B3/IL6/IL6RB/fine_screening/RoseTTAFold/protein_ligand/extracted \
  --protein-name IL6RB \
  --output-dir ~/Projects/B3/IL6/IL6RB/fine_screening/RoseTTAFold
```

### Step 2: Analyze Binding Sites

Run the analysis scripts to cluster binding sites and generate visualizations:

```bash
# AF3
python scripts/scoring/analyze_af3_binding_sites.py \
  --summary-csv ~/Projects/B3/IL6/IL6RB/fine_screening/AF3/summary.csv \
  --protein-name IL6RB \
  --output-dir ~/Projects/B3/IL6/IL6RB/fine_screening/AF3/binding_site_analysis \
  --cluster-eps 8.0

# Boltz2
python scripts/scoring/analyze_boltz2_binding_sites.py \
  --summary-csv ~/Projects/B3/IL6/IL6RB/fine_screening/Boltz2/summary.csv \
  --protein-name IL6RB \
  --output-dir ~/Projects/B3/IL6/IL6RB/fine_screening/Boltz2/binding_site_analysis \
  --cluster-eps 8.0

# RoseTTAFold
python scripts/scoring/analyze_rosettafold_binding_sites.py \
  --summary-csv ~/Projects/B3/IL6/IL6RB/fine_screening/RoseTTAFold/summary.csv \
  --protein-name IL6RB \
  --output-dir ~/Projects/B3/IL6/IL6RB/fine_screening/RoseTTAFold/binding_site_analysis \
  --cluster-eps 8.0
```

## Output Files

### Enhanced Summary CSV
The summary.csv files now include binding site information:
- Original columns (confidence metrics, SMILES, etc.)
- `binding_site_residues`: Residues within 5Å of ligand
- `binding_site_center`: Ligand center coordinates
- `num_binding_residues`: Number of binding site residues

### Analysis Output
Each analysis script generates:

1. **cluster_statistics.csv** - Per-cluster summary:
   - cluster_id: Cluster identifier (-1 for outliers)
   - num_compounds: Number of compounds in cluster
   - percentage: Percentage of total compounds
   - avg_ranking_score/avg_plddt: Average confidence metrics
   - avg_pae: Average PAE score
   - representative_residues: Top 10 most frequent residues

2. **binding_site_analysis.csv** - Full data with cluster assignments

3. **Visualizations**:
   - `binding_site_distribution.png` - Bar chart of cluster distribution
   - `binding_site_3d_scatter.png` - 3D scatter plot of binding site locations
   - `confidence_by_cluster.png` - Box plot of confidence by cluster
   - `residue_frequency_heatmap.png` - Heatmap of residue frequency per cluster

4. **PyMOL Script**:
   - `{protein_name}_binding_sites.pml` - Script to visualize binding sites
   - Colors binding site residues by cluster
   - Load with: `pymol {protein_name}_binding_sites.pml`

## Clustering Parameters

The analysis uses DBSCAN clustering with default parameters:
- `--cluster-eps 8.0`: Spatial tolerance in Angstroms (two sites within 8Å are neighbors)
- `--min-samples 3`: Minimum compounds to form a cluster

Adjust these parameters based on your protein size and binding site diversity.

## Implementation Details

### Binding Site Definition
- Distance cutoff: 5.0 Å (default)
- Measured from any protein atom to any ligand atom
- Only standard amino acid residues included (hetero residues excluded)

### Clustering Algorithm
- DBSCAN (Density-Based Spatial Clustering of Applications with Noise)
- Automatically discovers number of clusters
- Handles outliers (labeled as cluster -1)
- Clusters based on 3D coordinates of ligand center of mass

### Structure File Formats
- **AF3**: CIF files with protein (Chain A) and ligand (Chain Z)
- **Boltz2**: CIF files with protein (Chain A) and ligand (Chain Z, residue LIG1)
- **RoseTTAFold**: PDB files with protein (Chain A) and ligand (Chain B, residue LG1)

## Example Output

For a protein with ~3000 predictions, typical output might show:

```
cluster_id  num_compounds  percentage  avg_ranking_score  representative_residues
0           1847          61.9        0.78               A:CYS45,A:ARG67,A:GLU89,...
1           892           29.9        0.71               A:VAL23,A:CYS45,A:ARG67,...
2           156           5.2         0.65               A:SER156,A:ASP178,A:PHE201,...
-1          87            2.9         0.58               (outliers)
```

This reveals:
- Primary binding site (Cluster 0): 62% of compounds
- Secondary binding site (Cluster 1): 30% of compounds
- Tertiary binding site (Cluster 2): 5% of compounds
- Outliers: 3% with unique binding modes

## Troubleshooting

### Missing binding site information
- Check that structure files (.cif or .pdb) exist in the results folders
- Verify that ligand chains (Z for AF3/Boltz2, B for RoseTTAFold) are present
- Check for errors in the scoring script output

### No clusters found
- Try increasing `--cluster-eps` (e.g., 10.0 or 12.0)
- Try decreasing `--min-samples` (e.g., 2)
- Check that binding site centers are being extracted correctly

### Import errors
- Ensure all dependencies are installed: `pip install -r scripts/scoring/binding_site_requirements.txt`
- Check Python version (requires Python 3.7+)

## Notes

- The binding site extraction adds minimal overhead to the scoring scripts
- Analysis scripts can be run independently after scoring is complete
- Clustering parameters may need adjustment based on protein size and flexibility
- PyMOL scripts require a representative structure file to be loaded manually
