#!/usr/bin/env python3
"""Finalize peptide ranking from all AF3 results."""

import pandas as pd
import argparse
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="Aggregate and rank peptide validation results")
    parser.add_argument('--output-dir', required=True, help="Pipeline output directory")
    parser.add_argument('--protein-name', required=True, help="Protein name")

    args = parser.parse_args()

    results_dir = Path(args.output_dir) / "05_results"
    stream_dir = results_dir / "streaming_updates"
    mpnn_scores = Path(args.output_dir) / "02_proteinmpnn" / "scores.csv"

    # Aggregate all streaming batches
    all_dfs = []

    if stream_dir.exists():
        for batch_file in sorted(stream_dir.glob("batch_*.csv")):
            try:
                df = pd.read_csv(batch_file)
                all_dfs.append(df)
            except Exception as e:
                print(f"Warning: Could not read {batch_file}: {e}")

    if not all_dfs:
        print("Error: No streaming batch files found")
        return

    # Combine all results
    combined = pd.concat(all_dfs, ignore_index=True)

    # Join with ProteinMPNN scores if available
    if mpnn_scores.exists():
        try:
            mpnn_df = pd.read_csv(mpnn_scores)
            # Create matching key (assuming format matches)
            if 'backbone_id' in mpnn_df.columns:
                mpnn_df['join_key'] = mpnn_df['backbone_id'].astype(str) + '_' + mpnn_df['seq_idx'].astype(str)
                combined['join_key'] = combined['backbone_id'].astype(str) + '_seq' + combined.groupby('backbone_id').cumcount().add(1).astype(str)
                combined = combined.merge(mpnn_df[['join_key', 'score']], on='join_key', how='left')
                combined.rename(columns={'score': 'mpnn_score'}, inplace=True)
                combined.drop('join_key', axis=1, inplace=True)
        except Exception as e:
            print(f"Warning: Could not merge ProteinMPNN scores: {e}")

    # Sort by ranking score
    combined.sort_values('ranking_score', ascending=False, inplace=True)

    # Add rank column
    combined.insert(0, 'rank', range(1, len(combined) + 1))

    # Add AF3 confidence category
    def categorize_confidence(score):
        if score >= 0.7:
            return "high"
        elif score >= 0.5:
            return "medium"
        else:
            return "low"

    combined['af3_confidence'] = combined['ranking_score'].apply(categorize_confidence)

    # Save final results
    output_file = results_dir / "ranked_peptides.csv"
    combined.to_csv(output_file, index=False)

    print(f"\n{'='*60}")
    print(f"Peptide Design Results for {args.protein_name}")
    print(f"{'='*60}")
    print(f"Total peptides: {len(combined)}")
    print(f"High confidence: {(combined['af3_confidence'] == 'high').sum()}")
    print(f"Medium confidence: {(combined['af3_confidence'] == 'medium').sum()}")
    print(f"Low confidence: {(combined['af3_confidence'] == 'low').sum()}")
    print(f"\nTop 10 peptides:")
    print(combined[['rank', 'peptide_id', 'peptide_sequence', 'ranking_score', 'ipTM', 'pTM', 'pLDDT']].head(10).to_string(index=False))
    print(f"\nFull results saved to: {output_file}")
    print(f"{'='*60}\n")

if __name__ == "__main__":
    main()
