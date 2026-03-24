#!/usr/bin/env python

"""Utility for generating binding affinity predictions from trained HMSA checkpoints.

Example
-------
Assuming an input CSV named ``test_input.csv`` with columns ``SMILES`` and
``target`` (the latter is ignored during inference), run::

    python predict_binding.py \
        --test_path test_input.csv \
        --preds_path predictions.csv \
        --checkpoint_paths checkpoints/model.pt \
        --smiles_columns SMILES

This command generates ``predictions.csv`` with the predicted affinities and a
companion ``predictions_failed.csv`` that lists any rows filtered out during
preprocessing.
"""


from __future__ import annotations

import csv
import os
from typing import Iterable, List, Optional, Sequence, Tuple

import numpy as np
from tape import TAPETokenizer

from HMSA.args import PredictArgs
from HMSA.data import MoleculeDataLoader, MoleculeDatapoint, MoleculeDataset, get_data
from HMSA.train.predict import predict
from HMSA.utils import (
    load_args,
    load_checkpoint,
    load_scalers,
    makedirs,
    update_prediction_args,
)


class BindingPredictArgs(PredictArgs):
    """CLI arguments for the binding affinity prediction script."""

    failed_smiles_path: Optional[str] = None

    def process_args(self) -> None:
        super().process_args()

        if self.failed_smiles_path is None:
            base, ext = os.path.splitext(self.preds_path)
            if ext:
                self.failed_smiles_path = f"{base}_failed{ext}"
            else:
                self.failed_smiles_path = f"{self.preds_path}_failed.csv"


def _is_valid_datapoint(datapoint: MoleculeDatapoint) -> Tuple[bool, str]:
    """Determines whether the given datapoint can be processed by RDKit.

    Returns a tuple containing a boolean validity flag and a short failure reason.
    """

    try:
        if any(smiles == "" for smiles in datapoint.smiles):
            return False, "empty_smiles"

        mols = datapoint.mol
        for mol in mols:
            if isinstance(mol, tuple):
                reactant, product = mol
                if reactant is None or product is None:
                    return False, "rdkit_parse_failure"
                if reactant.GetNumHeavyAtoms() + product.GetNumHeavyAtoms() == 0:
                    return False, "no_heavy_atoms"
            else:
                if mol is None:
                    return False, "rdkit_parse_failure"
                if mol.GetNumHeavyAtoms() == 0:
                    return False, "no_heavy_atoms"
    except Exception as exc:  # pragma: no cover - defensive guard for RDKit edge cases.
        return False, f"exception:{exc.__class__.__name__}"

    return True, ""


def _write_failed_smiles(
    path: str,
    smiles_columns: Sequence[Optional[str]],
    failed: Sequence[Tuple[int, MoleculeDatapoint, str]],
) -> None:
    """Writes a CSV containing the SMILES strings that failed preprocessing."""

    makedirs(path, isfile=True)
    header = ["index"]
    header.extend(
        column if column is not None else f"smiles_{i}"
        for i, column in enumerate(smiles_columns)
    )
    header.append("reason")

    with open(path, "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(header)
        for index, datapoint, reason in failed:
            writer.writerow([index + 1, *datapoint.smiles, reason])


def _write_predictions(
    path: str,
    smiles_columns: Sequence[Optional[str]],
    task_names: Sequence[str],
    datapoints: Iterable[MoleculeDatapoint],
    predictions: np.ndarray,
    variances: Optional[np.ndarray] = None,
) -> None:
    """Writes a CSV containing SMILES strings alongside predicted affinities."""

    makedirs(path, isfile=True)

    smiles_header = [
        column if column is not None else f"smiles_{i}"
        for i, column in enumerate(smiles_columns)
    ]

    task_header = list(task_names)
    if not task_header:
        task_header = [f"affinity_{i + 1}" for i in range(predictions.shape[1])]

    header = smiles_header + task_header

    if variances is not None:
        variance_header = [f"{name}_variance" for name in task_header]
        header += variance_header

    with open(path, "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(header)

        for idx, datapoint in enumerate(datapoints):
            row = list(datapoint.smiles) + predictions[idx].tolist()
            if variances is not None:
                row += variances[idx].tolist()
            writer.writerow(row)


def _prepare_dataset(
    args: BindingPredictArgs,
) -> Tuple[MoleculeDataset, List[MoleculeDatapoint], List[Tuple[int, MoleculeDatapoint, str]]]:
    """Loads the raw data and separates valid and invalid datapoints."""

    raw_data = get_data(
        path=args.test_path,
        args=args,
        smiles_columns=args.smiles_columns,
        skip_invalid_smiles=False,
        store_row=True,
    )

    valid_datapoints: List[MoleculeDatapoint] = []
    failed_datapoints: List[Tuple[int, MoleculeDatapoint, str]] = []

    for index, datapoint in enumerate(raw_data):
        is_valid, reason = _is_valid_datapoint(datapoint)
        if is_valid:
            valid_datapoints.append(datapoint)
        else:
            failed_datapoints.append((index, datapoint, reason))

    dataset = MoleculeDataset(valid_datapoints)
    return dataset, valid_datapoints, failed_datapoints


def main() -> None:
    args = BindingPredictArgs().parse_args()

    if not args.checkpoint_paths:
        raise ValueError("At least one checkpoint must be provided for prediction.")

    primary_train_args = load_args(args.checkpoint_paths[0])
    update_prediction_args(args, primary_train_args)
    args.task_names = primary_train_args.task_names
    args.num_tasks = primary_train_args.num_tasks

    for extra_path in args.checkpoint_paths[1:]:
        extra_train_args = load_args(extra_path)
        if extra_train_args.task_names != primary_train_args.task_names:
            raise ValueError(
                "Checkpoint task definitions do not match; cannot aggregate predictions."
            )

    dataset, valid_datapoints, failed_datapoints = _prepare_dataset(args)

    _write_failed_smiles(args.failed_smiles_path, args.smiles_columns, failed_datapoints)

    if len(dataset) == 0:
        _write_predictions(
            path=args.preds_path,
            smiles_columns=args.smiles_columns,
            task_names=args.task_names or [],
            datapoints=valid_datapoints,
            predictions=np.zeros((0, args.num_tasks or 1)),
        )
        return

    data_loader = MoleculeDataLoader(
        dataset=dataset,
        batch_size=args.batch_size,
        num_workers=args.num_workers,
    )

    tokenizer = TAPETokenizer(vocab="unirep")

    sum_predictions = np.zeros((len(dataset), args.num_tasks), dtype=float)
    ensemble_predictions: Optional[List[np.ndarray]] = [] if args.ensemble_variance else None

    for checkpoint_path in args.checkpoint_paths:
        model = load_checkpoint(checkpoint_path, device=args.device)
        scaler, features_scaler, atom_descriptor_scaler, bond_feature_scaler = load_scalers(checkpoint_path)

        if features_scaler is not None:
            dataset.normalize_features(features_scaler)
        if atom_descriptor_scaler is not None:
            dataset.normalize_features(atom_descriptor_scaler, scale_atom_descriptors=True)
        if bond_feature_scaler is not None:
            dataset.normalize_features(bond_feature_scaler, scale_bond_features=True)

        model_predictions = np.array(
            predict(
                model=model,
                data_loader=data_loader,
                args=args,
                disable_progress_bar=True,
                scaler=scaler,
                tokenizer=tokenizer,
            )
        )

        sum_predictions += model_predictions

        if ensemble_predictions is not None:
            ensemble_predictions.append(model_predictions)

    mean_predictions = sum_predictions / len(args.checkpoint_paths)

    variance_matrix: Optional[np.ndarray] = None
    if ensemble_predictions is not None:
        stacked = np.stack(ensemble_predictions, axis=0)
        variance_matrix = stacked.var(axis=0)

    _write_predictions(
        path=args.preds_path,
        smiles_columns=args.smiles_columns,
        task_names=args.task_names or [],
        datapoints=valid_datapoints,
        predictions=mean_predictions,
        variances=variance_matrix,
    )


if __name__ == "__main__":
    main()
