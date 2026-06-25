"""Backfill binding_site_center / binding_site_residues / num_binding_residues
into an existing Boltz2 summary.csv.

Usage:
    python enrich_boltz2_summary.py --method-dir /path/to/fine_screening/Boltz2
"""
import argparse
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd
from tqdm import tqdm

sys.path.insert(0, os.path.dirname(__file__))
from boltz2_scores import extract_binding_site_residues, find_file


def resolve_cif(method_dir: str, folder: str) -> str | None:
    folder_path = os.path.join(method_dir, "output", folder)
    if not os.path.isdir(folder_path):
        return None
    return find_file(folder_path, "*_model_0.cif")


def process_row(args):
    method_dir, folder = args
    cif = resolve_cif(method_dir, folder)
    if cif is None:
        return folder, "", "", 0
    residues, center, count = extract_binding_site_residues(cif)
    return folder, residues, center, count


def main():
    parser = argparse.ArgumentParser(
        description="Enrich Boltz2 summary.csv with binding site coordinates",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--method-dir", required=True,
                        help="Path to fine_screening/Boltz2/ directory")
    parser.add_argument("--workers", type=int, default=os.cpu_count(),
                        help="Parallel worker processes")
    parser.add_argument("--force", action="store_true",
                        help="Re-process rows that already have binding_site_center")
    args = parser.parse_args()

    summary_path = os.path.join(args.method_dir, "summary.csv")
    if not os.path.exists(summary_path):
        print(f"Error: {summary_path} not found")
        sys.exit(1)

    df = pd.read_csv(summary_path)

    for col in ("binding_site_residues", "binding_site_center", "num_binding_residues"):
        if col not in df.columns:
            df[col] = "" if col != "num_binding_residues" else 0

    if args.force:
        to_process = df.index.tolist()
    else:
        to_process = df.index[df["binding_site_center"].isna() | (df["binding_site_center"] == "")].tolist()

    print(f"Rows to process: {len(to_process)} / {len(df)}")
    if not to_process:
        print("Nothing to do.")
        return

    tasks = [(args.method_dir, str(df.at[i, "folder"])) for i in to_process]

    results = {}
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(process_row, t): t[1] for t in tasks}
        for fut in tqdm(as_completed(futures), total=len(futures), desc="Enriching Boltz2"):
            folder, residues, center, count = fut.result()
            results[folder] = (residues, center, count)

    for i in to_process:
        folder = str(df.at[i, "folder"])
        if folder in results:
            residues, center, count = results[folder]
            df.at[i, "binding_site_residues"] = residues
            df.at[i, "binding_site_center"] = center
            df.at[i, "num_binding_residues"] = count

    tmp_path = summary_path + ".tmp"
    df.to_csv(tmp_path, index=False)
    os.replace(tmp_path, summary_path)
    enriched = df["binding_site_center"].notna() & (df["binding_site_center"] != "")
    print(f"Done. {enriched.sum()} / {len(df)} rows have binding_site_center.")


if __name__ == "__main__":
    main()
