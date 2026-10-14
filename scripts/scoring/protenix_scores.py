#!/usr/bin/env python3
"""Extract the highest ranked Protenix sample for each ligand."""

import argparse
from pathlib import Path

from ml_score_common import collect_method


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-folder", type=Path, required=True)
    parser.add_argument("--selected-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    path, counts = collect_method("Protenix", args.results_folder, args.selected_csv, args.output_dir)
    print(f"Wrote {path}: {counts}")


if __name__ == "__main__":
    main()
