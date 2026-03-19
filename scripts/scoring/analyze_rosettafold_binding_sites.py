#!/usr/bin/env python3
"""
Analyze binding sites from RoseTTAFold predictions.

This script reads the enhanced summary.csv file with binding site information
and performs clustering analysis to identify distinct binding modes.
"""

import os
import argparse
import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
from collections import Counter


def load_summary_with_binding_sites(summary_csv):
    """Load enhanced summary.csv with binding site information."""
    df = pd.read_csv(summary_csv)

    # Filter out rows with missing binding site information
    df = df[df['binding_site_center'].notna() & (df['binding_site_center'] != '')]

    # Parse binding site center coordinates
    def parse_coords(coord_str):
        try:
            x, y, z = map(float, coord_str.split(','))
            return np.array([x, y, z])
        except:
            return None

    df['coords'] = df['binding_site_center'].apply(parse_coords)
    df = df[df['coords'].notna()]

    return df


def cluster_binding_sites(df, eps=8.0, min_samples=3):
    """Cluster binding sites using DBSCAN based on spatial coordinates."""
    coords = np.vstack(df['coords'].values)

    # Perform DBSCAN clustering
    clustering = DBSCAN(eps=eps, min_samples=min_samples)
    labels = clustering.fit_predict(coords)

    df['cluster'] = labels
    return df


def calculate_cluster_statistics(df_with_clusters):
    """Calculate statistics for each cluster."""
    stats = []

    for cluster_id in sorted(df_with_clusters['cluster'].unique()):
        cluster_df = df_with_clusters[df_with_clusters['cluster'] == cluster_id]
        num_compounds = len(cluster_df)
        percentage = (num_compounds / len(df_with_clusters)) * 100

        # Calculate average confidence metrics (RoseTTAFold-specific)
        avg_plddt = cluster_df['mean_plddt'].mean()
        avg_pae = cluster_df['mean_pae'].mean()
        avg_pae_inter = cluster_df['pae_inter'].mean() if 'pae_inter' in cluster_df.columns else None

        # Find most common residues in this cluster
        all_residues = []
        for residues_str in cluster_df['binding_site_residues']:
            if pd.notna(residues_str) and residues_str:
                all_residues.extend(residues_str.split(','))

        residue_counts = Counter(all_residues)
        top_residues = [res for res, _ in residue_counts.most_common(10)]
        representative_residues = ','.join(top_residues)

        stats_dict = {
            'cluster_id': cluster_id,
            'num_compounds': num_compounds,
            'percentage': percentage,
            'avg_plddt': avg_plddt,
            'avg_pae': avg_pae,
            'representative_residues': representative_residues
        }
        if avg_pae_inter is not None:
            stats_dict['avg_pae_inter'] = avg_pae_inter
        stats.append(stats_dict)

    return pd.DataFrame(stats)


def visualize_binding_sites(df_with_clusters, output_dir, protein_name):
    """Generate visualizations of binding site clusters."""
    os.makedirs(output_dir, exist_ok=True)

    # 1. Binding site distribution (bar chart)
    cluster_counts = df_with_clusters['cluster'].value_counts().sort_index()
    percentages = (cluster_counts / len(df_with_clusters)) * 100

    plt.figure(figsize=(10, 6))
    bars = plt.bar(range(len(percentages)), percentages.values)
    plt.xlabel('Cluster ID', fontsize=12)
    plt.ylabel('Percentage of Compounds (%)', fontsize=12)
    plt.title(f'{protein_name} - Binding Site Distribution', fontsize=14)
    plt.xticks(range(len(percentages)), percentages.index)

    # Color outliers differently
    for i, cluster_id in enumerate(percentages.index):
        if cluster_id == -1:
            bars[i].set_color('gray')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'binding_site_distribution.png'), dpi=300)
    plt.close()

    # 2. 3D scatter plot
    coords = np.vstack(df_with_clusters['coords'].values)
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Color by cluster
    unique_clusters = sorted(df_with_clusters['cluster'].unique())
    colors = plt.cm.tab10(np.linspace(0, 1, len(unique_clusters)))

    for i, cluster_id in enumerate(unique_clusters):
        cluster_mask = df_with_clusters['cluster'] == cluster_id
        cluster_coords = coords[cluster_mask]
        cluster_scores = df_with_clusters[cluster_mask]['mean_plddt'].values

        label = f'Cluster {cluster_id}' if cluster_id != -1 else 'Outliers'
        ax.scatter(cluster_coords[:, 0], cluster_coords[:, 1], cluster_coords[:, 2],
                   c=[colors[i]], s=cluster_scores * 100, alpha=0.6, label=label)

    ax.set_xlabel('X (Å)', fontsize=12)
    ax.set_ylabel('Y (Å)', fontsize=12)
    ax.set_zlabel('Z (Å)', fontsize=12)
    ax.set_title(f'{protein_name} - Binding Site Spatial Distribution', fontsize=14)
    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'binding_site_3d_scatter.png'), dpi=300)
    plt.close()

    # 3. Confidence by cluster (box plot)
    plt.figure(figsize=(12, 6))
    cluster_data = []
    cluster_labels = []
    for cluster_id in sorted(df_with_clusters['cluster'].unique()):
        cluster_scores = df_with_clusters[df_with_clusters['cluster'] == cluster_id]['mean_plddt'].values
        cluster_data.append(cluster_scores)
        label = f'Cluster {cluster_id}' if cluster_id != -1 else 'Outliers'
        cluster_labels.append(label)

    plt.boxplot(cluster_data, labels=cluster_labels)
    plt.xlabel('Cluster', fontsize=12)
    plt.ylabel('Mean pLDDT', fontsize=12)
    plt.title(f'{protein_name} - Confidence by Cluster', fontsize=14)
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'confidence_by_cluster.png'), dpi=300)
    plt.close()

    # 4. Residue frequency heatmap
    # Get top 20 most frequent residues across all clusters
    all_residues = []
    for residues_str in df_with_clusters['binding_site_residues']:
        if pd.notna(residues_str) and residues_str:
            all_residues.extend(residues_str.split(','))

    residue_counts = Counter(all_residues)
    top_residues = [res for res, _ in residue_counts.most_common(20)]

    # Calculate frequency of each residue in each cluster
    heatmap_data = []
    cluster_ids = sorted([c for c in df_with_clusters['cluster'].unique() if c != -1])

    for residue in top_residues:
        row = []
        for cluster_id in cluster_ids:
            cluster_df = df_with_clusters[df_with_clusters['cluster'] == cluster_id]
            count = sum(residue in str(res_str) for res_str in cluster_df['binding_site_residues'])
            frequency = (count / len(cluster_df)) * 100 if len(cluster_df) > 0 else 0
            row.append(frequency)
        heatmap_data.append(row)

    if heatmap_data and cluster_ids:
        plt.figure(figsize=(10, 12))
        sns.heatmap(heatmap_data, xticklabels=[f'C{c}' for c in cluster_ids],
                    yticklabels=top_residues, cmap='YlOrRd', annot=False, cbar_kws={'label': 'Frequency (%)'})
        plt.xlabel('Cluster', fontsize=12)
        plt.ylabel('Residue', fontsize=12)
        plt.title(f'{protein_name} - Residue Frequency by Cluster', fontsize=14)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'residue_frequency_heatmap.png'), dpi=300)
        plt.close()


def generate_pymol_script(df_with_clusters, output_dir, protein_name, structure_path=None):
    """Generate PyMOL script to visualize binding sites."""
    script_path = os.path.join(output_dir, f'{protein_name}_binding_sites.pml')

    with open(script_path, 'w') as f:
        f.write(f"# PyMOL script for {protein_name} binding site visualization\n")
        f.write("# Generated by analyze_rosettafold_binding_sites.py\n\n")

        if structure_path:
            f.write(f"load {structure_path}\n")
        else:
            f.write("# Load your protein structure here\n")
            f.write("# load /path/to/structure.pdb\n\n")

        f.write("# Hide everything initially\n")
        f.write("hide everything\n")
        f.write("show cartoon, chain A\n")
        f.write("color gray80, chain A\n\n")

        # Get unique clusters (excluding outliers)
        clusters = sorted([c for c in df_with_clusters['cluster'].unique() if c != -1])
        colors = ['red', 'blue', 'green', 'yellow', 'orange', 'purple', 'cyan', 'magenta']

        for i, cluster_id in enumerate(clusters):
            cluster_df = df_with_clusters[df_with_clusters['cluster'] == cluster_id]
            color = colors[i % len(colors)]

            # Get all residues in this cluster
            all_residues = []
            for residues_str in cluster_df['binding_site_residues']:
                if pd.notna(residues_str) and residues_str:
                    all_residues.extend(residues_str.split(','))

            # Get unique residues
            unique_residues = list(set(all_residues))

            # Parse residue IDs
            resids = []
            for res in unique_residues:
                try:
                    # Format: A:CYS45 -> extract 45
                    resid = ''.join(filter(str.isdigit, res.split(':')[1]))
                    if resid:
                        resids.append(resid)
                except:
                    continue

            if resids:
                f.write(f"# Cluster {cluster_id} ({len(cluster_df)} compounds)\n")
                selection = f"cluster_{cluster_id}"
                resid_str = '+'.join(resids)
                f.write(f"select {selection}, chain A and resi {resid_str}\n")
                f.write(f"show sticks, {selection}\n")
                f.write(f"color {color}, {selection}\n\n")

        f.write("# Show ligand if present\n")
        f.write("show sticks, chain B\n")
        f.write("color white, chain B\n\n")

        f.write("# Set view\n")
        f.write("zoom chain A\n")
        f.write("bg_color white\n")

    print(f"PyMOL script saved to {script_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Analyze binding sites from RoseTTAFold predictions",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--summary-csv", type=str, required=True,
                        help="Path to enhanced summary.csv with binding site information")
    parser.add_argument("--protein-name", type=str, required=True,
                        help="Protein target name (e.g., IL6RB)")
    parser.add_argument("--output-dir", type=str, required=True,
                        help="Directory to save analysis results")
    parser.add_argument("--cluster-eps", type=float, default=8.0,
                        help="DBSCAN eps parameter (spatial tolerance in Angstroms)")
    parser.add_argument("--min-samples", type=int, default=3,
                        help="DBSCAN min_samples parameter")
    parser.add_argument("--structure-path", type=str, default=None,
                        help="Optional path to representative structure for PyMOL script")

    args = parser.parse_args()

    print(f"Loading summary from {args.summary_csv}")
    df = load_summary_with_binding_sites(args.summary_csv)
    print(f"Loaded {len(df)} predictions with binding site information")

    print(f"Clustering binding sites (eps={args.cluster_eps}, min_samples={args.min_samples})")
    df_clustered = cluster_binding_sites(df, eps=args.cluster_eps, min_samples=args.min_samples)

    print("Calculating cluster statistics")
    stats_df = calculate_cluster_statistics(df_clustered)

    # Save results
    os.makedirs(args.output_dir, exist_ok=True)
    stats_path = os.path.join(args.output_dir, 'cluster_statistics.csv')
    stats_df.to_csv(stats_path, index=False)
    print(f"Cluster statistics saved to {stats_path}")

    # Save full data with cluster assignments
    full_path = os.path.join(args.output_dir, 'binding_site_analysis.csv')
    df_clustered.drop(columns=['coords'], inplace=True)
    df_clustered.to_csv(full_path, index=False)
    print(f"Full analysis saved to {full_path}")

    print("Generating visualizations")
    visualize_binding_sites(df_clustered, args.output_dir, args.protein_name)

    print("Generating PyMOL script")
    generate_pymol_script(df_clustered, args.output_dir, args.protein_name, args.structure_path)

    print("\nAnalysis complete!")
    print(f"\nCluster Summary:")
    print(stats_df.to_string(index=False))


if __name__ == '__main__':
    main()

