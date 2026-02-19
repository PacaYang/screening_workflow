"""Utility script for running ColdDTA binding affinity predictions.

The script performs two stages:

1. Pre-process the input CSV by converting SMILES strings into graph
   representations and protein sequences into token indices. Any
   molecules that cannot be processed are collected and written to a
   failure report.
2. Load a trained ColdDTA checkpoint and generate affinity predictions
   for the successfully pre-processed molecules. The predictions are
   written to an output CSV.

Example
-------

```bash
python predict.py \
    --input data.csv \
    --checkpoint save/best_model.pt \
    --output predictions.csv \
    --failed-smiles failed.csv
```

The input CSV must contain the columns ``SMILES`` and
``sequence``. Extra columns are ignored.
"""

from __future__ import annotations

import argparse
import os
from typing import List, Sequence

import numpy as np
import pandas as pd
import torch
from rdkit import Chem
from torch_geometric import data as DATA
from torch_geometric.data import DataLoader
from tqdm import tqdm

from model import ColdDTA
from preprocessing import mol_to_graph_without_rm, seqs2int

TARGET_SEQUENCE_LEN = 1200


def _pad_or_truncate(sequence: Sequence[int], target_len: int) -> np.ndarray:
    """Ensure that the encoded sequence has a fixed length."""

    sequence = np.asarray(sequence, dtype=np.int64)
    if sequence.size < target_len:
        return np.pad(sequence, (0, target_len - sequence.size))
    return sequence[:target_len]


def preprocess_inputs(
    input_path: str,
    failed_output_path: str,
    target_len: int = TARGET_SEQUENCE_LEN,
) -> List[DATA.Data]:
    """Convert raw inputs into torch-geometric ``Data`` objects.

    Parameters
    ----------
    input_path:
        CSV file containing SMILES strings and protein target sequences.
    failed_output_path:
        Location where the SMILES that cannot be processed will be stored.
    target_len:
        Desired sequence length after padding/truncation.

    Returns
    -------
    List[Data]
        A list of pre-processed samples ready for inference.
    """

    df = pd.read_csv(input_path)

    data_list: List[DATA.Data] = []
    failures = []

    for _, row in tqdm(df.iterrows(), total=len(df), desc="Pre-processing"):
        smi = str(row.get("SMILES", "")).strip()
        sequence = str(row.get("sequence", "")).strip().upper()

        if not smi:
            failures.append({"SMILES": smi, "reason": "Missing SMILES"})
            continue
        if not sequence:
            failures.append({"SMILES": smi, "reason": "Missing target sequence"})
            continue

        try:
            mol = Chem.MolFromSmiles(smi)
        except Exception as exc:  # pragma: no cover - defensive branch
            failures.append({"SMILES": smi, "reason": f"RDKit error: {exc}"})
            continue

        if mol is None:
            failures.append({"SMILES": smi, "reason": "Invalid SMILES"})
            continue

        try:
            x, edge_index, edge_attr = mol_to_graph_without_rm(mol)
        except Exception as exc:
            failures.append({"SMILES": smi, "reason": f"Graph build error: {exc}"})
            continue

        try:
            target_tokens = seqs2int(sequence)
        except KeyError as exc:
            failures.append({"SMILES": smi, "reason": f"Unknown amino acid: {exc}"})
            continue

        target_tokens = _pad_or_truncate(target_tokens, target_len)

        data = DATA.Data(
            x=torch.FloatTensor(x),
            edge_index=edge_index,
            edge_attr=edge_attr,
            y=torch.zeros(1, dtype=torch.float32),
            target=torch.LongTensor([target_tokens]),
            smi=smi,
        )

        data_list.append(data)

    failure_df = pd.DataFrame(failures)
    if failure_df.empty:
        failure_df = pd.DataFrame(columns=["SMILES", "reason"])
    os.makedirs(os.path.dirname(failed_output_path) or ".", exist_ok=True)
    failure_df.to_csv(failed_output_path, index=False)

    return data_list


def predict_affinity(
    dataset: Sequence[DATA.Data],
    checkpoint_path: str,
    output_path: str,
    batch_size: int = 32,
) -> None:
    """Run model inference and save the predictions."""

    if not dataset:
        # Still create an empty prediction file to keep the workflow tidy.
        empty_df = pd.DataFrame(columns=["SMILES", "predicted_affinity"])
        os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
        empty_df.to_csv(output_path, index=False)
        return

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    model = ColdDTA().to(device)
    state_dict = torch.load(checkpoint_path, map_location=device)
    model.load_state_dict(state_dict)
    model.eval()

    loader = DataLoader(dataset, batch_size=batch_size, shuffle=False)

    smiles: List[str] = []
    predictions: List[float] = []

    for batch in tqdm(loader, desc="Predicting"):
        batch = batch.to(device)
        with torch.no_grad():
            output = model(batch).view(-1).cpu().numpy()

        batch_smiles = batch.smi
        if isinstance(batch_smiles, str):
            batch_smiles = [batch_smiles]
        else:
            batch_smiles = list(batch_smiles)

        smiles.extend(batch_smiles)
        predictions.extend(output.tolist())

    result_df = pd.DataFrame({
        "SMILES": smiles,
        "predicted_affinity": predictions,
    })
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    result_df.to_csv(output_path, index=False)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run ColdDTA predictions on new data.")
    parser.add_argument("--input", required=True, help="Path to the input CSV file.")
    parser.add_argument(
        "--checkpoint",
        required=True,
        help="Path to the trained ColdDTA checkpoint (state dict).",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Destination CSV for the predictions.",
    )
    parser.add_argument(
        "--failed-smiles",
        required=True,
        help="Destination CSV containing the SMILES strings that failed processing.",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=32,
        help="Batch size to use during inference.",
    )
    parser.add_argument(
        "--target-length",
        type=int,
        default=TARGET_SEQUENCE_LEN,
        help="Target sequence length after padding/truncation.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    dataset = preprocess_inputs(
        input_path=args.input,
        failed_output_path=args.failed_smiles,
        target_len=args.target_length,
    )

    predict_affinity(
        dataset=dataset,
        checkpoint_path=args.checkpoint,
        output_path=args.output,
        batch_size=args.batch_size,
    )


if __name__ == "__main__":
    main()

