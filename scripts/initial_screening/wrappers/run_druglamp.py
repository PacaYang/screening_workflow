#!/usr/bin/env python3
"""
DrugLAMP Wrapper Script

Converts workflow CSV format to DrugLAMP format, runs prediction, and converts output back.

Input format (from make_input_csv):
    SMILES, label, sequence

DrugLAMP input format:
    SMILES, Protein

DrugLAMP output format:
    SMILES, Protein, Prediction, Predicted_Label

Output format (for workflow):
    SMILES, sequence, druglamp_score
"""

import argparse
import pandas as pd
import subprocess
import sys
import tempfile
import os
from pathlib import Path


def convert_to_druglamp_format(input_csv, output_csv):
    """Convert workflow format to DrugLAMP format."""
    df = pd.read_csv(input_csv)

    # Check required columns
    if 'SMILES' not in df.columns or 'sequence' not in df.columns:
        raise ValueError(f"Input CSV must have 'SMILES' and 'sequence' columns. Found: {df.columns.tolist()}")

    # Create DrugLAMP format: SMILES, Protein
    druglamp_df = pd.DataFrame({
        'SMILES': df['SMILES'],
        'Protein': df['sequence']
    })

    druglamp_df.to_csv(output_csv, index=False)
    print(f"Converted {len(druglamp_df)} rows to DrugLAMP format")
    return len(druglamp_df)


def convert_druglamp_output(druglamp_csv, original_csv, output_csv):
    """Convert DrugLAMP output to workflow format."""
    druglamp_df = pd.read_csv(druglamp_csv)
    original_df = pd.read_csv(original_csv)

    # Check if DrugLAMP output has expected columns
    if 'Prediction' not in druglamp_df.columns:
        raise ValueError(f"DrugLAMP output missing 'Prediction' column. Found: {druglamp_df.columns.tolist()}")

    # Create output format: SMILES, sequence, druglamp_score
    output_df = pd.DataFrame({
        'SMILES': druglamp_df['SMILES'],
        'sequence': druglamp_df['Protein'],
        'druglamp_score': druglamp_df['Prediction']
    })

    output_df.to_csv(output_csv, index=False)
    print(f"Converted {len(output_df)} predictions to workflow format")
    return len(output_df)


def run_druglamp_inference(input_csv, output_csv, checkpoint, model, n_layer, device, batch_size, inference_script):
    """Run DrugLAMP inference."""
    # Create a simple temporary directory for inference
    # The fixed inference.py now properly handles paths internally
    temp_dir = tempfile.mkdtemp(prefix='druglamp_inference_')

    cmd = [
        'python', inference_script,
        '--checkpoint', checkpoint,
        '--input', input_csv,
        '--output', output_csv,
        '--model', model,
        '--n-layer', str(n_layer),
        '--device', device,
        '--batch-size', str(batch_size),
        '--temp-dir', temp_dir
    ]

    print(f"Running DrugLAMP inference: {' '.join(cmd)}")
    print(f"Working directory: {os.path.dirname(inference_script)}")
    print(f"Checkpoint file exists: {os.path.exists(checkpoint)}")
    print(f"Temp directory: {temp_dir}")

    try:
        # Don't capture output - let it stream to console for better debugging
        result = subprocess.run(cmd, check=True, text=True, cwd=os.path.dirname(inference_script))
        return True
    except subprocess.CalledProcessError as e:
        print(f"DrugLAMP inference failed with exit code {e.returncode}", file=sys.stderr)
        print(f"Command: {' '.join(cmd)}", file=sys.stderr)
        return False
    finally:
        # Clean up temporary directory
        if os.path.exists(temp_dir):
            import shutil
            shutil.rmtree(temp_dir, ignore_errors=True)


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
    parser = argparse.ArgumentParser(description="DrugLAMP wrapper for screening workflow")
    parser.add_argument('--input', required=True, help='Input CSV (workflow format)')
    parser.add_argument('--output', required=True, help='Output CSV (workflow format)')
    parser.add_argument('--failed-smiles', required=True, help='Failed SMILES CSV')
    parser.add_argument('--checkpoint', required=True, help='DrugLAMP checkpoint file')
    parser.add_argument('--inference-script', default='/home/yangl_pacagen_com/Applications/DrugLAMP/run_prediction_DrugLAMP.py',
                        help='Path to DrugLAMP inference.py')
    parser.add_argument('--model', default='DrugLAMP', help='Model architecture')
    parser.add_argument('--n-layer', type=int, default=30, help='ESM2 model size')
    parser.add_argument('--device', default='cuda', help='Device (cuda/cpu)')
    parser.add_argument('--batch-size', type=int, default=16, help='Batch size')

    args = parser.parse_args()

    # Create output directory
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.failed_smiles).parent.mkdir(parents=True, exist_ok=True)

    failed_smiles = []

    try:
        # Create temporary files for DrugLAMP format
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as tmp_input:
            tmp_input_path = tmp_input.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as tmp_output:
            tmp_output_path = tmp_output.name

        # Step 1: Convert input format
        print(f"Step 1: Converting input format...")
        n_input = convert_to_druglamp_format(args.input, tmp_input_path)

        # Step 2: Run DrugLAMP inference
        print(f"Step 2: Running DrugLAMP inference...")
        success = run_druglamp_inference(
            tmp_input_path, tmp_output_path,
            args.checkpoint, args.model, args.n_layer,
            args.device, args.batch_size, args.inference_script
        )

        if not success:
            print("DrugLAMP inference failed", file=sys.stderr)
            # Read input to get all SMILES as failed
            input_df = pd.read_csv(args.input)
            failed_smiles = input_df['SMILES'].tolist()
            write_failed_smiles(failed_smiles, args.failed_smiles)
            # Create empty output
            pd.DataFrame(columns=['SMILES', 'sequence', 'druglamp_score']).to_csv(args.output, index=False)
            sys.exit(1)

        # Step 3: Convert output format
        print(f"Step 3: Converting output format...")

        # Check if output file exists and has content
        if not os.path.exists(tmp_output_path):
            print(f"ERROR: DrugLAMP did not create output file: {tmp_output_path}", file=sys.stderr)
            # Read input to get all SMILES as failed
            input_df = pd.read_csv(args.input)
            failed_smiles = input_df['SMILES'].tolist()
            write_failed_smiles(failed_smiles, args.failed_smiles)
            # Create empty output
            pd.DataFrame(columns=['SMILES', 'sequence', 'druglamp_score']).to_csv(args.output, index=False)
            sys.exit(1)

        # Check if output file is empty
        output_size = os.path.getsize(tmp_output_path)
        print(f"DrugLAMP output file size: {output_size} bytes")
        if output_size == 0:
            print(f"ERROR: DrugLAMP output file is empty", file=sys.stderr)
            # Read input to get all SMILES as failed
            input_df = pd.read_csv(args.input)
            failed_smiles = input_df['SMILES'].tolist()
            write_failed_smiles(failed_smiles, args.failed_smiles)
            # Create empty output
            pd.DataFrame(columns=['SMILES', 'sequence', 'druglamp_score']).to_csv(args.output, index=False)
            sys.exit(1)

        n_output = convert_druglamp_output(tmp_output_path, args.input, args.output)

        # Check for failed SMILES (input count vs output count)
        if n_output < n_input:
            print(f"Warning: {n_input - n_output} SMILES failed during prediction")
            # Identify failed SMILES
            input_df = pd.read_csv(args.input)
            output_df = pd.read_csv(args.output)
            input_smiles = set(input_df['SMILES'])
            output_smiles = set(output_df['SMILES'])
            failed_smiles = list(input_smiles - output_smiles)

        write_failed_smiles(failed_smiles, args.failed_smiles)

        print(f"DrugLAMP wrapper completed successfully!")
        print(f"  Input: {n_input} compounds")
        print(f"  Output: {n_output} predictions")
        print(f"  Failed: {len(failed_smiles)} compounds")

    except Exception as e:
        print(f"Error in DrugLAMP wrapper: {e}", file=sys.stderr)
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
