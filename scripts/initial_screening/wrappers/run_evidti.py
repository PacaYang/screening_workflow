#!/usr/bin/env python3
"""
EviDTI Prediction Wrapper

Uses pre-extracted features to run EviDTI predictions.

Input format (from make_input_csv):
    SMILES, label, sequence

EviDTI input format:
    cid, uid, SMILES, seq, label

EviDTI output format:
    target_id, drug_id, prediction, probability, confidence, variance, true_label

Output format (for workflow):
    SMILES, sequence, evidti_score
"""

import argparse
import pandas as pd
import subprocess
import sys
import tempfile
import os
from pathlib import Path


def convert_to_evidti_format(input_csv, output_csv, protein_name):
    """Convert workflow format to EviDTI format."""
    df = pd.read_csv(input_csv)

    # Check required columns
    if 'SMILES' not in df.columns or 'sequence' not in df.columns:
        raise ValueError(f"Input CSV must have 'SMILES' and 'sequence' columns. Found: {df.columns.tolist()}")

    # Create EviDTI format: cid, uid, SMILES, seq, label
    evidti_df = pd.DataFrame({
        'cid': df['SMILES'],  # Use SMILES as compound ID
        'uid': [protein_name] * len(df),  # Protein name as UniProt ID
        'SMILES': df['SMILES'],
        'seq': df['sequence'],
        'label': [0] * len(df)  # Dummy label for prediction
    })

    evidti_df.to_csv(output_csv, index=False)
    print(f"Converted {len(evidti_df)} rows to EviDTI format")
    return len(evidti_df)


def convert_evidti_output(evidti_csv, original_csv, output_csv):
    """Convert EviDTI output to workflow format."""
    evidti_df = pd.read_csv(evidti_csv)
    original_df = pd.read_csv(original_csv)

    # Check if EviDTI output has expected columns
    if 'probability' not in evidti_df.columns:
        raise ValueError(f"EviDTI output missing 'probability' column. Found: {evidti_df.columns.tolist()}")

    # Create a mapping from drug_id (SMILES) to sequence
    smiles_to_seq = dict(zip(original_df['SMILES'], original_df['sequence']))

    # Create output format: SMILES, sequence, evidti_score
    # Use 'probability' as evidti_score (probability of positive class)
    output_df = pd.DataFrame({
        'SMILES': evidti_df['drug_id'],
        'sequence': evidti_df['drug_id'].map(smiles_to_seq),
        'evidti_score': evidti_df['probability']
    })

    # Remove rows where sequence mapping failed
    output_df = output_df.dropna(subset=['sequence'])

    output_df.to_csv(output_csv, index=False)
    print(f"Converted {len(output_df)} predictions to workflow format")
    return len(output_df)


def run_evidti_prediction(input_csv, output_csv, model_path, feature_dir, batch_size, predict_script):
    """Run EviDTI prediction."""
    cmd = [
        'python', predict_script,
        '--input', input_csv,
        '--output', output_csv,
        '--model_path', model_path,
        '--feature_dir', feature_dir,
        '--batch_size', str(batch_size),
        '--device', 'auto',
        '--no_labels'  # We're doing prediction only
    ]

    print(f"Running EviDTI prediction: {' '.join(cmd)}")

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(result.stdout)
        if result.stderr:
            print(f"EviDTI stderr: {result.stderr}", file=sys.stderr)
        return True
    except subprocess.CalledProcessError as e:
        print(f"EviDTI prediction failed: {e}", file=sys.stderr)
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
    parser = argparse.ArgumentParser(description="EviDTI wrapper for screening workflow")
    parser.add_argument('--input', required=True, help='Input CSV (workflow format)')
    parser.add_argument('--output', required=True, help='Output CSV (workflow format)')
    parser.add_argument('--failed-smiles', required=True, help='Failed SMILES CSV')
    parser.add_argument('--model-path', required=True, help='EviDTI model checkpoint')
    parser.add_argument('--feature-dir', required=True, help='Directory with pre-extracted features')
    parser.add_argument('--protein-name', required=True, help='Protein name/ID')
    parser.add_argument('--predict-script', default='/home/yangl_pacagen_com/Applications/EviDTI/predict_custom.py',
                        help='Path to EviDTI predict_custom.py')
    parser.add_argument('--batch-size', type=int, default=32, help='Batch size')

    args = parser.parse_args()

    # Create output directory
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.failed_smiles).parent.mkdir(parents=True, exist_ok=True)

    # Check if feature directory exists
    if not os.path.exists(args.feature_dir):
        print(f"Error: Feature directory not found: {args.feature_dir}", file=sys.stderr)
        sys.exit(1)

    # Check if required feature files exist
    required_features = ['protein_1d_feature.npy', 'drug_2d_feature.npy', 'drug_3d_feature.npy']
    for feature_file in required_features:
        feature_path = os.path.join(args.feature_dir, feature_file)
        if not os.path.exists(feature_path):
            print(f"Error: Required feature file not found: {feature_path}", file=sys.stderr)
            sys.exit(1)

    failed_smiles = []

    try:
        # Create temporary files for EviDTI format
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as tmp_input:
            tmp_input_path = tmp_input.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as tmp_output:
            tmp_output_path = tmp_output.name

        # Step 1: Convert input format
        print(f"Step 1: Converting input format...")
        n_input = convert_to_evidti_format(args.input, tmp_input_path, args.protein_name)

        # Step 2: Run EviDTI prediction
        print(f"Step 2: Running EviDTI prediction...")
        success = run_evidti_prediction(
            tmp_input_path, tmp_output_path,
            args.model_path, args.feature_dir, args.batch_size, args.predict_script
        )

        if not success:
            print("EviDTI prediction failed", file=sys.stderr)
            # All SMILES failed
            input_df = pd.read_csv(args.input)
            failed_smiles = input_df['SMILES'].tolist()
            write_failed_smiles(failed_smiles, args.failed_smiles)
            # Create empty output
            pd.DataFrame(columns=['SMILES', 'sequence', 'evidti_score']).to_csv(args.output, index=False)
            sys.exit(1)

        # Step 3: Convert output format
        print(f"Step 3: Converting output format...")
        n_output = convert_evidti_output(tmp_output_path, args.input, args.output)

        # Check for failed SMILES
        if n_output < n_input:
            print(f"Warning: {n_input - n_output} SMILES failed during prediction")
            # Identify failed SMILES
            input_df = pd.read_csv(args.input)
            output_df = pd.read_csv(args.output)
            input_smiles = set(input_df['SMILES'])
            output_smiles = set(output_df['SMILES'])
            failed_smiles = list(input_smiles - output_smiles)

        write_failed_smiles(failed_smiles, args.failed_smiles)

        print(f"EviDTI wrapper completed successfully!")
        print(f"  Input: {n_input} compounds")
        print(f"  Output: {n_output} predictions")
        print(f"  Failed: {len(failed_smiles)} compounds")

    except Exception as e:
        print(f"Error in EviDTI wrapper: {e}", file=sys.stderr)
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
