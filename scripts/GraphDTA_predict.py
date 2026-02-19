"""Utility script for running GraphDTA predictions on new data.

This script expects an input CSV file containing at least two columns:
```
compound_iso_smiles,target_sequence
```
It will attempt to convert every SMILES string into the graph
representation required by the GraphDTA models.  Any SMILES strings that
cannot be parsed by RDKit are written to a failure file.  Predictions for
successfully parsed entries are written to the requested output file.
"""
import argparse
import csv
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from rdkit import Chem
from rdkit.Chem import Mol
from torch_geometric import data as DATA
from torch_geometric.data import DataLoader

from models.gat import GATNet
from models.gat_gcn import GAT_GCN
from models.gcn import GCNNet
from models.ginconv import GINConvNet

SEQ_VOCAB = "ABCDEFGHIKLMNOPQRSTUVWXYZ"
SEQ_DICT = {v: (i + 1) for i, v in enumerate(SEQ_VOCAB)}
MAX_SEQ_LEN = 1000


class InvalidSequenceError(ValueError):
    """Raised when a protein sequence contains no valid amino acid codes."""


def one_of_k_encoding(x: int, allowable_set: Sequence[int]) -> List[bool]:
    if x not in allowable_set:
        raise ValueError(f"Input {x} not in allowable set {allowable_set}")
    return [x == s for s in allowable_set]


def one_of_k_encoding_unk(x: int, allowable_set: Sequence[int]) -> List[bool]:
    """Map inputs not present in the allowable set to the last element."""

    if x not in allowable_set:
        x = allowable_set[-1]
    return [x == s for s in allowable_set]


def atom_features(atom: Chem.Atom) -> np.ndarray:
    features = np.array(
        one_of_k_encoding_unk(
            atom.GetSymbol(),
            [
                "C",
                "N",
                "O",
                "S",
                "F",
                "Si",
                "P",
                "Cl",
                "Br",
                "Mg",
                "Na",
                "Ca",
                "Fe",
                "As",
                "Al",
                "I",
                "B",
                "V",
                "K",
                "Tl",
                "Yb",
                "Sb",
                "Sn",
                "Ag",
                "Pd",
                "Co",
                "Se",
                "Ti",
                "Zn",
                "H",
                "Li",
                "Ge",
                "Cu",
                "Au",
                "Ni",
                "Cd",
                "In",
                "Mn",
                "Zr",
                "Cr",
                "Pt",
                "Hg",
                "Pb",
                "Unknown",
            ],
        )
        + one_of_k_encoding(atom.GetDegree(), list(range(11)))
        + one_of_k_encoding_unk(atom.GetTotalNumHs(), list(range(11)))
        + one_of_k_encoding_unk(atom.GetImplicitValence(), list(range(11)))
        + [atom.GetIsAromatic()]
    )
    denom = features.sum()
    if denom == 0:
        return features
    return features / denom


def mol_to_graph(mol: Mol) -> Tuple[int, List[List[float]], List[List[int]]]:
    """Convert an RDKit Mol into graph data."""

    c_size = mol.GetNumAtoms()
    features = [atom_features(atom) for atom in mol.GetAtoms()]

    edges = []
    for bond in mol.GetBonds():
        edges.append([bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()])

    graph = []
    for u, v in edges:
        graph.append([u, v])
        graph.append([v, u])

    if not graph:
        graph = [[0, 0]]

    return c_size, features, graph


def seq_to_tensor(seq: str) -> torch.LongTensor:
    seq_indices: List[int] = []
    for ch in seq[:MAX_SEQ_LEN]:
        seq_indices.append(SEQ_DICT.get(ch, 0))

    if not seq_indices:
        raise InvalidSequenceError("Empty protein sequence after preprocessing")

    if len(seq_indices) < MAX_SEQ_LEN:
        seq_indices.extend([0] * (MAX_SEQ_LEN - len(seq_indices)))

    return torch.LongTensor([seq_indices])


def build_data_object(
    smile_graph: Tuple[int, Sequence[Sequence[float]], Sequence[Sequence[int]]],
    seq_tensor: torch.LongTensor,
) -> DATA.Data:
    c_size, features, edge_index = smile_graph
    data = DATA.Data(
        x=torch.tensor(features, dtype=torch.float32),
        edge_index=torch.tensor(edge_index, dtype=torch.long).t().contiguous(),
        y=torch.zeros(1, dtype=torch.float32),
    )
    data.target = seq_tensor.clone().long()
    data.__setitem__("c_size", torch.tensor([c_size], dtype=torch.long))
    return data


def load_model(model_type: str, weights_path: Path, device: torch.device) -> nn.Module:
    model_map = {
        "gin": GINConvNet,
        "gat": GATNet,
        "gat_gcn": GAT_GCN,
        "gcn": GCNNet,
    }
    try:
        model_class = model_map[model_type.lower()]
    except KeyError as exc:
        raise ValueError(
            f"Unknown model type '{model_type}'. Choose from {sorted(model_map)}"
        ) from exc

    model = model_class().to(device)
    state_dict = torch.load(weights_path, map_location=device)
    model.load_state_dict(state_dict)
    model.eval()
    return model


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run GraphDTA predictions")
    parser.add_argument("input", type=Path, help="Path to the input CSV file")
    parser.add_argument(
        "--model",
        type=Path,
        required=True,
        help="Path to the trained model weights (.model file)",
    )
    parser.add_argument(
        "--model-type",
        default="gat",
        choices=["gin", "gat", "gat_gcn", "gcn"],
        help="Model architecture used during training",
    )
    parser.add_argument(
        "--batch-size", type=int, default=512, help="Mini-batch size for inference"
    )
    parser.add_argument(
        "--device",
        default="cuda" if torch.cuda.is_available() else "cpu",
        help="Computation device (e.g. 'cpu', 'cuda', 'cuda:0')",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("predictions.csv"),
        help="File to write predictions to",
    )
    parser.add_argument(
        "--failed-smiles",
        type=Path,
        default=Path("failed_smiles.txt"),
        help="File to record SMILES strings that cannot be processed",
    )
    parser.add_argument(
        "--smiles-column",
        default="compound_iso_smiles",
        help="Column name containing SMILES strings",
    )
    parser.add_argument(
        "--sequence-column",
        default="sequence",
        help="Column name containing protein sequences",
    )
    return parser.parse_args()


def write_failed_smiles(
    failed_smiles: Iterable[str],
    path: Path,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["smiles"])
        for smiles in failed_smiles:
            writer.writerow([smiles])


def main() -> None:
    args = parse_args()

    df = pd.read_csv(args.input)

    failed_smiles: List[str] = []
    data_objects: List[DATA.Data] = []
    metadata: List[Tuple[str, str]] = []

    for idx, row in df.iterrows():
        smiles = str(row[args.smiles_column]).strip()
        sequence = str(row[args.sequence_column]).strip()

        if not smiles:
            failed_smiles.append(smiles)
            continue

        mol = Chem.MolFromSmiles(smiles)
        if mol is None or mol.GetNumAtoms() == 0:
            failed_smiles.append(smiles)
            continue

        try:
            smile_graph = mol_to_graph(mol)
            seq_tensor = seq_to_tensor(sequence)
        except InvalidSequenceError:
            failed_smiles.append(smiles)
            continue

        data = build_data_object(smile_graph, seq_tensor)
        data_objects.append(data)
        metadata.append((smiles, sequence))

    write_failed_smiles(failed_smiles, args.failed_smiles)

    if not data_objects:
        print("No valid records were found in the input file. Nothing to predict.")
        return

    device = torch.device(args.device)
    model = load_model(args.model_type, args.model, device)

    loader = DataLoader(data_objects, batch_size=args.batch_size, shuffle=False)

    predictions: List[float] = []
    with torch.no_grad():
        for batch in loader:
            batch = batch.to(device)
            output = model(batch)
            predictions.extend(output.view(-1).cpu().tolist())

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            ["SMILES", "sequence", "predicted_affinity"]
        )
        for (smiles, sequence), pred in zip(metadata, predictions):
            writer.writerow([smiles, sequence, pred])

    print(
        f"Wrote {len(predictions)} predictions to {args.output} and "
        f"{len(failed_smiles)} failing SMILES to {args.failed_smiles}."
    )


if __name__ == "__main__":
    main()
