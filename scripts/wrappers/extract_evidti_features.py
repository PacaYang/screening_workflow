#!/usr/bin/env python3
"""
EviDTI Feature Extraction Wrapper

Orchestrates the 3-step feature extraction process for EviDTI:
1. Extract protein 1D features (ProtT5 embeddings)
2. Extract drug 2D features (Molecular graph embeddings)
3. Extract drug 3D features (3D conformer features)

Then creates aligned feature arrays for prediction.

This script is designed to be run once per protein, with features cached
and reused across all input chunks.
"""

import argparse
import pandas as pd
import subprocess
import sys
import os
from pathlib import Path
import glob
import tempfile
import shutil


def create_modified_script(original_script, output_dir, protein_sequences_csv, drug_info_csv, smiles_dir):
    """Create a temporary modified copy of an EviDTI script with correct paths."""
    with open(original_script, 'r') as f:
        content = f.read()

    script_name = os.path.basename(original_script)

    # Modify paths based on which script this is
    if 'extract_p_emb.py' in script_name:
        # Replace hardcoded dataset path
        content = content.replace(
            "load_dataset('dataset/bisai/target_info.csv')",
            f"load_dataset('{protein_sequences_csv}')"
        )
        # Replace hardcoded save path
        content = content.replace(
            "'dataset/bisai/protein_emb.npy'",
            f"'{os.path.join(output_dir, 'protein_1d_feature.npy')}'"
        )
        # Replace Windows path for ProtT5 tokenizer and model
        content = content.replace(
            "'D:\\zyp\\hub\\models--Rostlab--prot_t5_xl_uniref50\\snapshots\\973be27c52ee6474de9c945952a8008aeb2a1a73'",
            "'Rostlab/prot_t5_xl_uniref50'"
        )
        # Also handle double backslashes (in case of different escaping)
        content = content.replace(
            "'D:\\\\zyp\\\\hub\\\\models--Rostlab--prot_t5_xl_uniref50\\\\snapshots\\\\973be27c52ee6474de9c945952a8008aeb2a1a73'",
            "'Rostlab/prot_t5_xl_uniref50'"
        )

    elif 'extract_drug_2d_emb.py' in script_name:
        # Replace hardcoded dataset paths
        content = content.replace(
            "pd.read_csv('dataset/bisai/drug_info.csv')['ID']",
            f"pd.read_csv('{drug_info_csv}')['ID']"
        )
        content = content.replace(
            "Graph_Classification_Dataset('dataset/bisai/drug_info.csv'",
            f"Graph_Classification_Dataset('{drug_info_csv}'"
        )
        # Replace hardcoded save path
        content = content.replace(
            "'dataset/bisai/drug_2d_feature_raw.npy'",
            f"'{os.path.join(output_dir, 'drug_2d_feature.npy')}'"
        )

    elif 'extract_drug_3d_emb.py' in script_name:
        # Replace hardcoded save path
        content = content.replace(
            "'dataset/bisai/drug_3d_feature_raw.npy'",
            f"'{os.path.join(output_dir, 'drug_3d_feature.npy')}'"
        )

    # Create temporary file
    temp_script = tempfile.NamedTemporaryFile(mode='w', suffix='.py', delete=False, dir=os.path.dirname(original_script))
    temp_script.write(content)
    temp_script.close()

    return temp_script.name


def combine_input_files(input_dir):
    """Combine all input_*.csv files from the checkpoint."""
    csv_files = sorted(glob.glob(os.path.join(input_dir, "input_*.csv")))

    if not csv_files:
        raise ValueError(f"No input_*.csv files found in {input_dir}")

    print(f"Found {len(csv_files)} input files to combine")

    dfs = []
    for csv_file in csv_files:
        df = pd.read_csv(csv_file)
        dfs.append(df)

    combined_df = pd.concat(dfs, ignore_index=True)
    print(f"Combined {len(combined_df)} total rows")

    return combined_df


def prepare_evidti_input(combined_df, output_dir, protein_name):
    """Prepare input CSV in EviDTI format."""
    # EviDTI expects: cid, uid, SMILES, seq, label (optional)
    evidti_df = pd.DataFrame({
        'cid': combined_df['SMILES'],  # Use SMILES as compound ID
        'uid': [protein_name] * len(combined_df),  # Protein name as UniProt ID
        'SMILES': combined_df['SMILES'],
        'seq': combined_df['sequence'],
        'label': [0] * len(combined_df)  # Dummy label for prediction
    })

    # Remove duplicates
    evidti_df = evidti_df.drop_duplicates(subset=['SMILES', 'seq'])

    input_csv = os.path.join(output_dir, 'input_validated.csv')
    evidti_df.to_csv(input_csv, index=False)
    print(f"Prepared EviDTI input: {len(evidti_df)} unique drug-target pairs")

    return input_csv


def run_feature_preparation(input_csv, output_dir, prepare_script):
    """Run prepare_features_custom.py to set up feature extraction."""
    cmd = [
        'python', prepare_script,
        '--input', input_csv,
        '--output_dir', output_dir
    ]

    print(f"Running feature preparation: {' '.join(cmd)}")

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True, cwd=os.path.dirname(prepare_script))
        print(result.stdout)
        if result.stderr:
            print(f"Preparation stderr: {result.stderr}", file=sys.stderr)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Feature preparation failed: {e}", file=sys.stderr)
        print(f"stdout: {e.stdout}", file=sys.stderr)
        print(f"stderr: {e.stderr}", file=sys.stderr)
        return False


def run_protein_feature_extraction(output_dir, extract_p_script):
    """Extract protein 1D features using ProtT5."""
    print("Extracting protein 1D features (ProtT5)...")

    # Create paths for the files that prepare_features_custom.py created
    protein_sequences_csv = os.path.join(output_dir, 'protein_sequences.csv')
    drug_info_csv = os.path.join(output_dir, 'drug_info.csv')
    smiles_dir = os.path.join(output_dir, 'smiles')

    # Create modified script with correct paths
    temp_script = create_modified_script(extract_p_script, output_dir, protein_sequences_csv, drug_info_csv, smiles_dir)

    try:
        cmd = ['python', temp_script]
        print(f"Running: {' '.join(cmd)}")

        # Don't capture output - let it stream for better debugging
        result = subprocess.run(cmd, check=True, text=True,
                                cwd=os.path.dirname(extract_p_script))

        # Check if output file was created
        protein_feature_file = os.path.join(output_dir, 'protein_1d_feature.npy')
        if os.path.exists(protein_feature_file):
            print(f"✓ Protein features extracted: {protein_feature_file}")
            return True
        else:
            print(f"✗ Protein feature file not found: {protein_feature_file}", file=sys.stderr)
            return False

    except subprocess.CalledProcessError as e:
        print(f"Protein feature extraction failed with exit code {e.returncode}", file=sys.stderr)
        return False
    finally:
        # Clean up temporary script
        if os.path.exists(temp_script):
            os.unlink(temp_script)


def run_drug_2d_feature_extraction(output_dir, extract_2d_script):
    """Extract drug 2D features using Molecular BERT."""
    print("Extracting drug 2D features (Molecular BERT)...")

    # Create paths for the files that prepare_features_custom.py created
    protein_sequences_csv = os.path.join(output_dir, 'protein_sequences.csv')
    drug_info_csv = os.path.join(output_dir, 'drug_info.csv')
    smiles_dir = os.path.join(output_dir, 'smiles')

    # Create modified script with correct paths
    temp_script = create_modified_script(extract_2d_script, output_dir, protein_sequences_csv, drug_info_csv, smiles_dir)

    try:
        cmd = ['python', temp_script]
        print(f"Running: {' '.join(cmd)}")

        # Don't capture output - let it stream for better debugging
        result = subprocess.run(cmd, check=True, text=True,
                                cwd=os.path.dirname(extract_2d_script))

        # Check if output file was created
        drug_2d_feature_file = os.path.join(output_dir, 'drug_2d_feature.npy')
        if os.path.exists(drug_2d_feature_file):
            print(f"✓ Drug 2D features extracted: {drug_2d_feature_file}")
            return True
        else:
            print(f"✗ Drug 2D feature file not found: {drug_2d_feature_file}", file=sys.stderr)
            return False

    except subprocess.CalledProcessError as e:
        print(f"Drug 2D feature extraction failed with exit code {e.returncode}", file=sys.stderr)
        return False
    finally:
        # Clean up temporary script
        if os.path.exists(temp_script):
            os.unlink(temp_script)


def run_drug_3d_feature_extraction(output_dir, extract_3d_script):
    """Extract drug 3D features using GeoGNN."""
    print("Extracting drug 3D features (GeoGNN)...")

    # The 3D extraction script may need the smiles directory
    smiles_dir = os.path.join(output_dir, 'smiles')

    # Create paths for the files that prepare_features_custom.py created
    protein_sequences_csv = os.path.join(output_dir, 'protein_sequences.csv')
    drug_info_csv = os.path.join(output_dir, 'drug_info.csv')

    # Create modified script with correct paths
    temp_script = create_modified_script(extract_3d_script, output_dir, protein_sequences_csv, drug_info_csv, smiles_dir)

    try:
        cmd = ['python', temp_script, '--data_path', smiles_dir]
        print(f"Running: {' '.join(cmd)}")

        # Don't capture output - let it stream for better debugging
        result = subprocess.run(cmd, check=True, text=True,
                                cwd=os.path.dirname(extract_3d_script))

        # Check if output file was created
        drug_3d_feature_file = os.path.join(output_dir, 'drug_3d_feature.npy')
        if os.path.exists(drug_3d_feature_file):
            print(f"✓ Drug 3D features extracted: {drug_3d_feature_file}")
            return True
        else:
            print(f"✗ Drug 3D feature file not found: {drug_3d_feature_file}", file=sys.stderr)
            return False

    except subprocess.CalledProcessError as e:
        print(f"Drug 3D feature extraction failed with exit code {e.returncode}", file=sys.stderr)
        return False
    finally:
        # Clean up temporary script
        if os.path.exists(temp_script):
            os.unlink(temp_script)


def create_final_features(input_csv, output_dir, prepare_script):
    """Create aligned feature arrays using prepare_features_custom.py --create_final."""
    cmd = [
        'python', prepare_script,
        '--input', input_csv,
        '--output_dir', output_dir,
        '--create_final'
    ]

    print(f"Creating final aligned features: {' '.join(cmd)}")

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True,
                                cwd=os.path.dirname(prepare_script))
        print(result.stdout)
        if result.stderr:
            print(f"Final feature creation stderr: {result.stderr}", file=sys.stderr)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Final feature creation failed: {e}", file=sys.stderr)
        print(f"stdout: {e.stdout}", file=sys.stderr)
        print(f"stderr: {e.stderr}", file=sys.stderr)
        return False


def main():
    parser = argparse.ArgumentParser(description="EviDTI feature extraction wrapper")
    parser.add_argument('--input-dir', required=True, help='Directory with input_*.csv files')
    parser.add_argument('--output-dir', required=True, help='Output directory for features')
    parser.add_argument('--protein-name', required=True, help='Protein name/ID')
    parser.add_argument('--extract-p-script', default='/home/yangl_pacagen_com/Applications/EviDTI/extract_p_emb.py')
    parser.add_argument('--extract-2d-script', default='/home/yangl_pacagen_com/Applications/EviDTI/extract_drug_2d_emb.py')
    parser.add_argument('--extract-3d-script', default='/home/yangl_pacagen_com/Applications/EviDTI/extract_drug_3d_emb.py')
    parser.add_argument('--prepare-script', default='/home/yangl_pacagen_com/Applications/EviDTI/prepare_features_custom.py')

    args = parser.parse_args()

    # Create output directory
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    try:
        # Step 1: Combine all input files
        print("=" * 60)
        print("Step 1: Combining input files...")
        print("=" * 60)
        combined_df = combine_input_files(args.input_dir)

        # Step 2: Prepare EviDTI input format
        print("\n" + "=" * 60)
        print("Step 2: Preparing EviDTI input format...")
        print("=" * 60)
        input_csv = prepare_evidti_input(combined_df, args.output_dir, args.protein_name)

        # Step 3: Run feature preparation
        print("\n" + "=" * 60)
        print("Step 3: Running feature preparation...")
        print("=" * 60)
        if not run_feature_preparation(input_csv, args.output_dir, args.prepare_script):
            print("Feature preparation failed", file=sys.stderr)
            sys.exit(1)

        # Step 4: Extract protein features
        print("\n" + "=" * 60)
        print("Step 4: Extracting protein 1D features...")
        print("=" * 60)
        if not run_protein_feature_extraction(args.output_dir, args.extract_p_script):
            print("Protein feature extraction failed", file=sys.stderr)
            sys.exit(1)

        # Step 5: Extract drug 2D features
        print("\n" + "=" * 60)
        print("Step 5: Extracting drug 2D features...")
        print("=" * 60)
        if not run_drug_2d_feature_extraction(args.output_dir, args.extract_2d_script):
            print("Drug 2D feature extraction failed", file=sys.stderr)
            sys.exit(1)

        # Step 6: Extract drug 3D features
        print("\n" + "=" * 60)
        print("Step 6: Extracting drug 3D features...")
        print("=" * 60)
        if not run_drug_3d_feature_extraction(args.output_dir, args.extract_3d_script):
            print("Drug 3D feature extraction failed", file=sys.stderr)
            sys.exit(1)

        # Step 7: Create final aligned features
        print("\n" + "=" * 60)
        print("Step 7: Creating final aligned features...")
        print("=" * 60)
        if not create_final_features(input_csv, args.output_dir, args.prepare_script):
            print("Final feature creation failed", file=sys.stderr)
            sys.exit(1)

        # Create done token
        done_token = os.path.join(args.output_dir, 'feature_extraction.done')
        Path(done_token).touch()

        print("\n" + "=" * 60)
        print("✓ Feature extraction completed successfully!")
        print("=" * 60)
        print(f"Features saved to: {args.output_dir}")
        print(f"  - protein_1d_feature.npy")
        print(f"  - drug_2d_feature.npy")
        print(f"  - drug_3d_feature.npy")

    except Exception as e:
        print(f"Error in feature extraction: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
