#!/usr/bin/env python3
"""
ConPLex Wrapper Script

Converts workflow CSV format to ConPLex TSV format, runs prediction, and converts output back.

Input format (from make_input_csv):
    SMILES, label, sequence

ConPLex input format (TSV, no header):
    proteinID, moleculeID, proteinSequence, moleculeSMILES

ConPLex output format (TSV, no header):
    moleculeID, proteinID, prediction_score

Output format (for workflow):
    SMILES, sequence, conplex_score
"""

import argparse
import pandas as pd
import subprocess
import sys
import tempfile
import os
from pathlib import Path


def convert_to_conplex_format(input_csv, output_tsv, protein_name):
    """Convert workflow format to ConPLex TSV format."""
    df = pd.read_csv(input_csv)

    # Check required columns
    if 'SMILES' not in df.columns or 'sequence' not in df.columns:
        raise ValueError(f"Input CSV must have 'SMILES' and 'sequence' columns. Found: {df.columns.tolist()}")

    # Create ConPLex format: proteinID, moleculeID, proteinSequence, moleculeSMILES
    # Use SMILES as moleculeID for easy mapping back
    conplex_df = pd.DataFrame({
        'proteinID': [protein_name] * len(df),
        'moleculeID': df['SMILES'],  # Use SMILES as ID
        'proteinSequence': df['sequence'],
        'moleculeSMILES': df['SMILES']
    })

    # Write TSV without header
    conplex_df.to_csv(output_tsv, sep='\t', header=False, index=False)
    print(f"Converted {len(conplex_df)} rows to ConPLex TSV format")
    return df


def convert_conplex_output(conplex_tsv, original_df, output_csv):
    """Convert ConPLex TSV output to workflow CSV format."""
    # Read ConPLex output (TSV, no header)
    conplex_df = pd.read_csv(conplex_tsv, sep='\t', header=None,
                              names=['moleculeID', 'proteinID', 'prediction_score'])

    # Map moleculeID (SMILES) back to original data
    # Create a mapping from SMILES to sequence
    smiles_to_seq = dict(zip(original_df['SMILES'], original_df['sequence']))

    # Create output format: SMILES, sequence, conplex_score
    output_df = pd.DataFrame({
        'SMILES': conplex_df['moleculeID'],
        'sequence': conplex_df['moleculeID'].map(smiles_to_seq),
        'conplex_score': conplex_df['prediction_score']
    })

    # Remove rows where sequence mapping failed
    output_df = output_df.dropna(subset=['sequence'])

    output_df.to_csv(output_csv, index=False)
    print(f"Converted {len(output_df)} predictions to workflow format")
    return len(output_df)


def run_conplex_prediction(input_tsv, output_tsv, model_path, device, batch_size, predict_script):
    """Run ConPLex prediction."""
    cmd = [
        'python', predict_script,
        '--input', input_tsv,
        '--output', output_tsv,
        '--model-path', model_path,
        '--device', device,
        '--batch-size', str(batch_size)
    ]

    print(f"Running ConPLex prediction: {' '.join(cmd)}")

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(result.stdout)
        if result.stderr:
            print(f"ConPLex stderr: {result.stderr}", file=sys.stderr)
        return True
    except subprocess.CalledProcessError as e:
        print(f"ConPLex prediction failed: {e}", file=sys.stderr)
        print(f"stdout: {e.stdout}", file=sys.stderr)
        print(f"stderr: {e.stderr}", file=sys.stderr)
        return False


def write_failed_smiles(failed_smiles, output_path):
    """Write failed SMILES to CSV."""
    if failed_smiles:
        failed_df = pd.DataFrame({'SMILES': failed_smiles})
        failed_df.to_csv(output_path, index=False)
        print(f"Wrote {len(failed_smiles)} failed SMILES to {output_path}")
    else:
        # Create empty file
        Path(output_path).touch()
        print(f"No failed SMILES, created empty file at {output_path}")


def main():
    parser = argparse.ArgumentParser(description="ConPLex wrapper for screening workflow")
    parser.add_argument('--input', required=True, help='Input CSV (workflow format)')
    parser.add_argument('--output', required=True, help='Output CSV (workflow format)')
    parser.add_argument('--failed-smiles', required=True, help='Failed SMILES CSV')
    parser.add_argument('--model-path', required=True, help='ConPLex model file')
    parser.add_argument('--protein-name', required=True, help='Protein name/ID')
    parser.add_argument('--predict-script', default='/home/yangl_pacagen_com/screening_workflow/scripts/run_prediction_ConPLex.py',
                        help='Path to ConPLex run_prediction_jobs.py')
    parser.add_argument('--device', default='0', help='GPU device ID or cpu')
    parser.add_argument('--batch-size', type=int, default=128, help='Batch size')

    args = parser.parse_args()

    # Create output directory
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.failed_smiles).parent.mkdir(parents=True, exist_ok=True)

    failed_smiles = []

    try:
        # Create temporary files for ConPLex format
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as tmp_input:
            tmp_input_path = tmp_input.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as tmp_output:
            tmp_output_path = tmp_output.name

        # Step 1: Convert input format
        print(f"Step 1: Converting input format...")
        original_df = convert_to_conplex_format(args.input, tmp_input_path, args.protein_name)
        n_input = len(original_df)

        # Step 2: Run ConPLex prediction
        print(f"Step 2: Running ConPLex prediction...")
        success = run_conplex_prediction(
            tmp_input_path, tmp_output_path,
            args.model_path, args.device, args.batch_size, args.predict_script
        )

        if not success:
            print("ConPLex prediction failed", file=sys.stderr)
            # All SMILES failed
            failed_smiles = original_df['SMILES'].tolist()
            write_failed_smiles(failed_smiles, args.failed_smiles)
            # Create empty output
            pd.DataFrame(columns=['SMILES', 'sequence', 'conplex_score']).to_csv(args.output, index=False)
            sys.exit(1)

        # Step 3: Convert output format
        print(f"Step 3: Converting output format...")
        n_output = convert_conplex_output(tmp_output_path, original_df, args.output)

        # Check for failed SMILES
        if n_output < n_input:
            print(f"Warning: {n_input - n_output} SMILES failed during prediction")
            # Identify failed SMILES
            output_df = pd.read_csv(args.output)
            input_smiles = set(original_df['SMILES'])
            output_smiles = set(output_df['SMILES'])
            failed_smiles = list(input_smiles - output_smiles)

        write_failed_smiles(failed_smiles, args.failed_smiles)

        print(f"ConPLex wrapper completed successfully!")
        print(f"  Input: {n_input} compounds")
        print(f"  Output: {n_output} predictions")
        print(f"  Failed: {len(failed_smiles)} compounds")

    except Exception as e:
        print(f"Error in ConPLex wrapper: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)

    finally:
        # Clean up temporary files
        if 'tmp_input_path' in locals() and os.path.exists(tmp_input_path):
            os.unlink(tmp_input_path)
        if 'tmp_output_path' in locals() and os.path.exists(tmp_output_path):
            os.unlink(tmp_output_path)


if __name__ == '__main__':
    main()
