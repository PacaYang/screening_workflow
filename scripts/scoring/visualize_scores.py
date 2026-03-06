import argparse
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Draw
from io import BytesIO
from PIL import Image


def load_af3(target_dir):
    """Load AF3 summary.csv and return (df, score_col, label, title)."""
    path = os.path.join(target_dir, "fine_screening", "AF3", "summary.csv")
    if not os.path.exists(path):
        return None, None, None, None
    df = pd.read_csv(path)
    score_col = "chain_pair_pae_min"
    # Filter out -1 values (failed predictions)
    df = df[df[score_col] != -1].copy()
    df = df.sort_values(score_col, ascending=True).reset_index(drop=True)
    return df, score_col, "PAE (A)", "AF3"


def load_boltz2(target_dir):
    """Load Boltz2 summary.csv and return (df, score_col, label, title)."""
    path = os.path.join(target_dir, "fine_screening", "Boltz2", "summary.csv")
    if not os.path.exists(path):
        return None, None, None, None
    df = pd.read_csv(path)
    score_col = "affinity"
    df = df.sort_values(score_col, ascending=True).reset_index(drop=True)
    return df, score_col, "Affinity (kcal/mol)", "Boltz2"


def load_vina(target_dir):
    """Load Vina results.csv, take best affinity per SMILES, return (df, score_col, label, title)."""
    path = os.path.join(target_dir, "fine_screening", "Vina", "results.csv")
    if not os.path.exists(path):
        return None, None, None, None
    df = pd.read_csv(path)
    score_col = "affinity"
    # Take best (most negative) affinity per SMILES across all boxes
    df = df.loc[df.groupby("SMILES")[score_col].idxmin()].copy()
    df = df.sort_values(score_col, ascending=True).reset_index(drop=True)
    return df, score_col, "Affinity (kcal/mol)", "Vina"


def smiles_to_image(smiles, size=(250, 250)):
    """Convert a SMILES string to a PIL Image."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        # Return a blank image with error text
        img = Image.new('RGB', size, 'white')
        return img
    img = Draw.MolToImage(mol, size=size)
    return img


def plot_distributions(datasets, target_name, output_dir):
    """Create score distribution histograms for each algorithm."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle(f"{target_name} - Score Distributions", fontsize=16, fontweight='bold')

    for idx, (df, score_col, ylabel, algo_name) in enumerate(datasets):
        ax = axes[idx]
        if df is None or df.empty:
            ax.text(0.5, 0.5, f"{algo_name}\nNo data available",
                    ha='center', va='center', transform=ax.transAxes, fontsize=14)
            ax.set_title(algo_name)
            continue

        scores = df[score_col].dropna()
        mean_val = scores.mean()
        median_val = scores.median()

        ax.hist(scores, bins=50, color='steelblue', edgecolor='black', alpha=0.7)
        ax.axvline(mean_val, color='red', linestyle='--', linewidth=1.5, label=f'Mean: {mean_val:.2f}')
        ax.axvline(median_val, color='orange', linestyle='-', linewidth=1.5, label=f'Median: {median_val:.2f}')
        ax.set_title(algo_name, fontsize=14)
        ax.set_xlabel(ylabel, fontsize=12)
        ax.set_ylabel("Count", fontsize=12)
        ax.legend(fontsize=10)

    plt.tight_layout()
    out_path = os.path.join(output_dir, f"{target_name}_score_distributions.png")
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {out_path}")


def plot_top5(datasets, target_name, output_dir):
    """Create a grid of top-5 compounds with 2D structures for each algorithm."""
    n_rows = 5
    n_cols = 3
    cell_w, cell_h = 3.0, 3.5
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(cell_w * n_cols, cell_h * n_rows))
    fig.suptitle(f"{target_name} - Top 5 Compounds", fontsize=16, fontweight='bold', y=1.02)

    for col_idx, (df, score_col, ylabel, algo_name) in enumerate(datasets):
        for row_idx in range(n_rows):
            ax = axes[row_idx, col_idx]
            ax.set_xticks([])
            ax.set_yticks([])

            if row_idx == 0:
                ax.set_title(algo_name, fontsize=13, fontweight='bold')

            if df is None or row_idx >= len(df):
                ax.text(0.5, 0.5, "N/A", ha='center', va='center',
                        transform=ax.transAxes, fontsize=12, color='gray')
                continue

            row = df.iloc[row_idx]
            smiles = row["SMILES"]
            score = row[score_col]

            img = smiles_to_image(smiles, size=(300, 300))
            ax.imshow(img)
            ax.set_xlabel(f"{ylabel}: {score:.3f}", fontsize=9)

    plt.tight_layout()
    out_path = os.path.join(output_dir, f"{target_name}_top5_compounds.png")
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {out_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Visualize fine screening scores for a target",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--target-dir", type=str, required=True,
                        help="Path to the target directory (e.g., ~/Projects/B3/IL6/IL6RA)")
    parser.add_argument("--target-name", type=str, required=True,
                        help="Target name for plot titles (e.g., IL6RA)")
    args = parser.parse_args()

    target_dir = os.path.expanduser(args.target_dir)

    # Load data from all three algorithms
    af3_data = load_af3(target_dir)
    boltz2_data = load_boltz2(target_dir)
    vina_data = load_vina(target_dir)

    datasets = [af3_data, boltz2_data, vina_data]

    # Generate plots
    plot_distributions(datasets, args.target_name, target_dir)
    plot_top5(datasets, args.target_name, target_dir)


if __name__ == "__main__":
    main()
