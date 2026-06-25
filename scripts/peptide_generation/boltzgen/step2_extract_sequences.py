#!/usr/bin/env python3
"""Extract designed peptide sequences from BoltzGen output.

BoltzGen's filter step writes:
  <output>/diverse_<budget>/final_designs_metrics_<budget>.csv

The CSV contains 'designed_chain_sequence' (full chain) and 'designed_sequence'
(designed residues only). We use 'designed_chain_sequence' as the peptide
sequence for AF3 validation, falling back to 'designed_sequence' if absent.

Output: one .fa file per design in <output_dir>/seqs/, compatible with
step3_prepare_af3.py --no-skip-native.
"""

import argparse
import csv
import sys
from pathlib import Path


def find_metrics_csv(boltzgen_output: Path) -> Path:
    """Locate the final_designs_metrics CSV in the boltzgen output tree."""
    # Primary location: diverse_<N>/final_designs_metrics_<N>.csv
    candidates = sorted(boltzgen_output.glob("diverse_*/final_designs_metrics_*.csv"))
    if candidates:
        return candidates[-1]  # take the largest budget if multiple

    # Fallback: any final_designs_metrics CSV anywhere in the output
    candidates = sorted(boltzgen_output.rglob("final_designs_metrics_*.csv"))
    if candidates:
        return candidates[-1]

    # Last resort: all_designs_metrics.csv (unfiltered)
    fallback = boltzgen_output / "all_designs_metrics.csv"
    if fallback.exists():
        print(
            "Warning: final_designs_metrics not found, falling back to all_designs_metrics.csv",
            file=sys.stderr,
        )
        return fallback

    raise FileNotFoundError(
        f"No BoltzGen metrics CSV found under {boltzgen_output}. "
        "Ensure BoltzGen completed the filter step."
    )


def extract_sequences(boltzgen_output: Path, output_dir: Path) -> int:
    metrics_csv = find_metrics_csv(boltzgen_output)
    print(f"Reading designs from: {metrics_csv}")

    seqs_dir = output_dir / "seqs"
    seqs_dir.mkdir(parents=True, exist_ok=True)

    count = 0
    with open(metrics_csv, newline="") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames or []

        seq_col = (
            "designed_chain_sequence"
            if "designed_chain_sequence" in fieldnames
            else "designed_sequence"
        )
        id_col = "id" if "id" in fieldnames else None

        for i, row in enumerate(reader):
            seq = row.get(seq_col, "").strip()
            if not seq:
                continue

            design_id = row.get(id_col, f"design_{i}") if id_col else f"design_{i}"
            # Sanitise id for use in filenames
            safe_id = str(design_id).replace("/", "_").replace(" ", "_")

            fa_path = seqs_dir / f"{safe_id}.fa"
            fa_path.write_text(f">{safe_id}\n{seq}\n")
            count += 1

    return count


def main():
    parser = argparse.ArgumentParser(
        description="Extract BoltzGen designed sequences into per-design .fa files"
    )
    parser.add_argument(
        "--boltzgen-output",
        required=True,
        help="Path to boltzgen output directory (01_boltzgen/output)",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Output directory (sequences written to <output-dir>/seqs/)",
    )
    args = parser.parse_args()

    boltzgen_output = Path(args.boltzgen_output)
    output_dir = Path(args.output_dir)

    if not boltzgen_output.exists():
        print(f"Error: BoltzGen output directory not found: {boltzgen_output}", file=sys.stderr)
        sys.exit(1)

    n = extract_sequences(boltzgen_output, output_dir)

    if n == 0:
        print("Error: No sequences extracted. Check BoltzGen output.", file=sys.stderr)
        sys.exit(1)

    print(f"Extracted {n} sequences to {output_dir}/seqs/")


if __name__ == "__main__":
    main()
