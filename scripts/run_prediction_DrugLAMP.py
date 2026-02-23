"""
DrugLAMP Inference Script
Usage:
    python inference.py --checkpoint <path_to_checkpoint> --input <input_csv> --output <output_csv>

Input CSV format:
    SMILES,Protein
    CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O,MKKFFDSRREQGGSGLGSGSSGGGGSTSGLGSGYIGRVFGIGRQQVTVDEVLAEGGFAIVFLVRTSNGMKCALKRMFVNNEHDLQVCKREIQIMRDLSGHKNIVGYIDSSINNVSSGDVWEVLILMDFCRGGQVVNLMNQRLQTGFTENEVLQIFCDTCEAVARLHQCKTPIIHRDLKVENILLHDRGHYVLCDFGSATNKFQNPQTEGVNAVEDEIKKYTTLSYRAPEMVNLYSGKIITTKADIWALGCLLYKLCYFTLPFGESQVAICDGNFTIPDNSRYSQDMHCLIRYMLEPDPDKRPDIYQVSYFSFKLLKKECPIPNVQNSPIPAKLPEPVKASEAAAKKTQPKARLTDPIPTTETSIAPRQRPKAGQTQPNQAQKPPQAPPPQGQKVTPPPPQAPSPNQAPPPPPQTQPPPPQPSPPASPIPPQAPQQPVEHQPPPPPRPQKPQAPPPPGQVTPIPPPLPQAPAQPEVQPPPPPTPIPPQVPQAPPAQPPPPSQAPPQQSQPPQPPQAQHQQPQPQHQQPPPQHPPPPQQPPQHQHQQPPPQQHQHQHQPPPQQPPPHQHPPQHQQPQQHPPQPHHHQQPPQQQPPQQHPPPHPPQPQPHQPPQQPQQSPPQQHQQQQPQQPQPQPPPPQQHHQQHQPQPPQPQPHHQHQQQQQHQPPPHHPQHQPPPPQPHQQPQPPQQQPQPQQQPPQHQPPPQQHQPPHPPPQSPPPPKPQLPPPHHQQHQQPPPPPPPPPQPQPQPPPPHHHQPQHPPQPPQQQQQHHPPPPQQPQPQPQPQQQHPPPQP
"""

import os
import sys
import torch
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from tqdm import tqdm

# Add DrugLAMP to path
script_dir = Path(__file__).parent
sys.path.insert(0, str(script_dir))

import esm.pretrained as esp
from configs import get_cfg_defaults
from model import MInterface
from handler import MultiModalityDataset
from utils import set_seed, multimodality_collate_func
from torch.utils.data import DataLoader

# ESM2 model options
n_layer2esp_fns = {
    48: esp.esm2_t48_15B_UR50D,
    36: esp.esm2_t36_3B_UR50D,
    33: esp.esm2_t33_650M_UR50D,
    30: esp.esm2_t30_150M_UR50D,
    12: esp.esm2_t12_35M_UR50D
}


class DrugLAMPInference:
    def __init__(self, checkpoint_path, model_name='DrugLAMP', n_layer=30, device='cuda'):
        """
        Initialize DrugLAMP for inference

        Args:
            checkpoint_path: Path to trained model checkpoint (.ckpt file)
            model_name: Model architecture name (default: 'DrugLAMP')
            n_layer: ESM2 model size (12, 30, 33, 36, or 48)
            device: Device to run inference on ('cuda' or 'cpu')
        """
        self.checkpoint_path = checkpoint_path
        self.model_name = model_name
        self.n_layer = n_layer
        self.device = torch.device(device if torch.cuda.is_available() else 'cpu')

        print(f"Using device: {self.device}")

        # Load configuration
        model_cfg = script_dir / 'configs' / f'{model_name}.yaml'
        self.cfg = get_cfg_defaults()
        if model_cfg.exists():
            self.cfg.merge_from_file(str(model_cfg))

        set_seed(self.cfg.SOLVER.SEED)

        # Load ESM2 model
        print(f"Loading ESM2 model (n_layer={n_layer})...")
        self.esp_fn = n_layer2esp_fns[n_layer]

        # Model will be loaded later with dataset parameters
        self.model = None

    def load_model(self, n_drug_feature, n_prot_feature):
        """Load the trained model"""
        print("Loading DrugLAMP model...")

        # Create model interface
        model_interface = MInterface(self.model_name, self.cfg)
        self.model = model_interface.load_model(
            n_drug_feature=n_drug_feature,
            n_prot_feature=n_prot_feature
        )

        # Load checkpoint weights
        if os.path.exists(self.checkpoint_path):
            print(f"Loading checkpoint from {self.checkpoint_path}")
            checkpoint = torch.load(self.checkpoint_path, map_location=self.device)

            # Extract state dict (handle PyTorch Lightning checkpoint format)
            if 'state_dict' in checkpoint:
                state_dict = checkpoint['state_dict']
                # Remove 'exp_model.' prefix if present
                state_dict = {k.replace('exp_model.', ''): v for k, v in state_dict.items()}
            else:
                state_dict = checkpoint

            self.model.load_state_dict(state_dict, strict=False)
            print("Checkpoint loaded successfully!")
        else:
            raise FileNotFoundError(f"Checkpoint not found: {self.checkpoint_path}")

        self.model.to(self.device)
        self.model.eval()

    def prepare_data(self, input_csv, temp_dir='temp_inference'):
        """
        Prepare input data for inference

        Args:
            input_csv: Path to CSV file with columns: SMILES, Protein
            temp_dir: Temporary directory for processed data
        """
        # Create temporary directory structure
        os.makedirs(temp_dir, exist_ok=True)

        # Read input data
        df = pd.read_csv(input_csv)
        required_cols = ['SMILES', 'Protein']

        if not all(col in df.columns for col in required_cols):
            raise ValueError(f"Input CSV must contain columns: {required_cols}")

        # Add a dummy label column (required by dataset class)
        # The dataset expects a column named 'Y' for the label
        if 'Y' not in df.columns:
            df['Y'] = 0

        # Save as full.csv and test.csv in temp directory
        df.to_csv(os.path.join(temp_dir, 'full.csv'), index=False)
        df.to_csv(os.path.join(temp_dir, 'test.csv'), index=False)

        # Create dataset
        print("Preparing dataset...")
        max_drug_atoms = self.cfg.DRUG.MAX_NODES
        dataset = MultiModalityDataset(
            temp_dir,
            'test.csv',
            self.esp_fn,
            self.n_layer,
            self.device,
            gen_embed=True,
            max_drug_atoms=max_drug_atoms
        )

        return dataset

    def predict(self, input_csv, output_csv=None, batch_size=1, temp_dir='temp_inference'):
        """
        Run inference on input data

        Args:
            input_csv: Path to input CSV file
            output_csv: Path to save predictions (optional)
            batch_size: Batch size for inference
            temp_dir: Temporary directory for processing
        """
        # Prepare dataset
        dataset = self.prepare_data(input_csv, temp_dir)

        # Load model with dataset parameters
        if self.model is None:
            self.load_model(dataset.n_drug_feature, dataset.n_prot_feature)

        # Create dataloader
        dataloader = DataLoader(
            dataset,
            batch_size=batch_size,
            shuffle=False,
            num_workers=0,
            collate_fn=multimodality_collate_func
        )

        # Run inference
        print("Running inference...")
        predictions = []
        smiles_list = []
        protein_list = []

        with torch.no_grad():
            for batch in tqdm(dataloader, desc="Processing batches"):
                feat_d, feat_p, labels, llm_d, llm_p, meta = batch

                # Move to device
                feat_d = feat_d.to(self.device)
                feat_p = feat_p.to(self.device)
                llm_d = llm_d.to(self.device)
                llm_p = llm_p.to(self.device)

                # Forward pass
                _, _, _, _, score = self.model(feat_d, feat_p, llm_d, llm_p, mode='train')

                # Apply sigmoid for binary classification
                probs = torch.sigmoid(score).cpu().numpy()
                predictions.extend(probs.flatten().tolist())

                # Collect SMILES and Protein from meta
                for m in meta:
                    smiles_list.append(m['Drug'])
                    protein_list.append(m['Prot'])

        # Create dataframe with only valid predictions (failed SMILES are excluded)
        df = pd.DataFrame({
            'SMILES': smiles_list,
            'Protein': protein_list,
            'Prediction': predictions
        })
        df['Predicted_Label'] = (df['Prediction'] > 0.5).astype(int)

        print(f"Generated {len(df)} predictions")
        if len(dataset.failed_drug_ords) > 0:
            print(f"Note: {len(dataset.failed_drug_ords)} SMILES were skipped due to processing errors")

        # Save results
        if output_csv:
            df.to_csv(output_csv, index=False)
            print(f"Results saved to {output_csv}")

        return df


def main():
    parser = argparse.ArgumentParser(description="DrugLAMP Inference Script")
    parser.add_argument('--checkpoint', required=True, type=str,
                        help='Path to trained model checkpoint (.ckpt file)')
    parser.add_argument('--input', required=True, type=str,
                        help='Path to input CSV file (columns: SMILES, Protein)')
    parser.add_argument('--output', type=str, default='predictions.csv',
                        help='Path to output CSV file with predictions')
    parser.add_argument('--model', type=str, default='DrugLAMP',
                        choices=['DrugLAMP', 'DrugLAMP2C2P', 'DrugLAMPwoLLM'],
                        help='Model architecture name')
    parser.add_argument('--n-layer', type=int, default=30,
                        choices=[12, 30, 33, 36, 48],
                        help='ESM2 model size (number of layers)')
    parser.add_argument('--device', type=str, default='cuda',
                        choices=['cuda', 'cpu'],
                        help='Device to run inference on')
    parser.add_argument('--batch-size', type=int, default=1,
                        help='Batch size for inference')
    parser.add_argument('--temp-dir', type=str, default='temp_inference',
                        help='Temporary directory for processing')

    args = parser.parse_args()

    # Validate inputs
    if not os.path.exists(args.checkpoint):
        print(f"Error: Checkpoint file not found: {args.checkpoint}")
        sys.exit(1)

    if not os.path.exists(args.input):
        print(f"Error: Input file not found: {args.input}")
        sys.exit(1)

    # Run inference
    inferencer = DrugLAMPInference(
        checkpoint_path=args.checkpoint,
        model_name=args.model,
        n_layer=args.n_layer,
        device=args.device
    )

    results = inferencer.predict(
        input_csv=args.input,
        output_csv=args.output,
        batch_size=args.batch_size,
        temp_dir=args.temp_dir
    )

    print("\nInference Summary:")
    print(f"Total predictions: {len(results)}")
    print(f"Positive predictions: {results['Predicted_Label'].sum()}")
    print(f"Negative predictions: {len(results) - results['Predicted_Label'].sum()}")
    print(f"\nPrediction statistics:")
    print(results['Prediction'].describe())


if __name__ == '__main__':
    main()
