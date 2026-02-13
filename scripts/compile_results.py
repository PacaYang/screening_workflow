import os
import argparse
import pandas as pd

BASE = "/shared/B3/IL4_13_all"

# Per-algorithm config: (subfolder to summary CSV, ranking column, ascending?)
ALGO_CONFIG = {
    "AF3": {
        "csv_path": "fine_screening/AF3/summary.csv",
        "rank_col": "chain_pair_pae_min",
        "ascending": True,
        "filter_val": -1,  # filter out rows where rank_col == this value
    },
    "Boltz2": {
        "csv_path": "fine_screening/Boltz2/summary.csv",
        "rank_col": "affinity",
        "ascending": True,
    },
    "Vina": {
        "csv_path": "fine_screening/Vina/results.csv",
        "rank_col": "affinity",
        "ascending": True,
    },
    "RoseTTAFold": {
        "csv_path": "fine_screening/RoseTTAFold/summary.csv",
        "rank_col": "pae_inter",
        "ascending": True,
    },
}

# Which algorithms are expected per target
TARGET_ALGOS = {
    "IL13":    ["AF3", "Boltz2", "Vina"],
    "IL13RA1": ["AF3", "Boltz2"],
    "IL4":     ["AF3", "Boltz2", "Vina"],
    "IL4RA":   ["AF3", "Boltz2", "Vina", "RoseTTAFold"],
}

TOP_N = 150


def load_and_rank(target, algo, top_n=TOP_N):
    cfg = ALGO_CONFIG[algo]
    csv_path = os.path.join(BASE, target, cfg["csv_path"])

    if not os.path.exists(csv_path):
        print(f"  WARNING: {csv_path} not found, skipping {algo} for {target}")
        return None

    df = pd.read_csv(csv_path)

    # Ensure SMILES column exists
    if "SMILES" not in df.columns:
        print(f"  WARNING: No SMILES column in {csv_path}, skipping")
        return None

    rank_col = cfg["rank_col"]
    if rank_col not in df.columns:
        print(f"  WARNING: Column '{rank_col}' not in {csv_path}, skipping")
        return None

    # Convert rank column to numeric
    df[rank_col] = pd.to_numeric(df[rank_col], errors="coerce")

    # Filter out invalid rows
    filter_val = cfg.get("filter_val")
    if filter_val is not None:
        df = df[df[rank_col] != filter_val]

    df = df.dropna(subset=[rank_col])

    # Sort and take top N
    df = df.sort_values(rank_col, ascending=cfg["ascending"]).head(top_n).copy()
    df["algorithm"] = algo
    df["rank_in_algorithm"] = range(1, len(df) + 1)

    return df


def compile_target(target, output_dir):
    algos = TARGET_ALGOS.get(target, [])
    frames = []

    for algo in algos:
        print(f"  Processing {algo}...")
        df = load_and_rank(target, algo)
        if df is not None and len(df) > 0:
            frames.append(df)
            print(f"    -> {len(df)} rows")

    if not frames:
        print(f"  No data for {target}")
        return

    combined = pd.concat(frames, ignore_index=True)

    # Reorder columns: SMILES, algorithm, rank_in_algorithm, then everything else
    priority_cols = ["SMILES", "algorithm", "rank_in_algorithm"]
    other_cols = [c for c in combined.columns if c not in priority_cols]
    combined = combined[priority_cols + other_cols]

    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, f"{target}_compiled.csv")
    combined.to_csv(out_path, index=False)
    print(f"  Wrote {len(combined)} rows to {out_path}")


def main():
    parser = argparse.ArgumentParser(description="Compile fine screening results")
    parser.add_argument("--output-dir", type=str,
                        default=os.path.join(BASE, "compiled_results"),
                        help="Directory for compiled output CSVs")
    parser.add_argument("--targets", nargs="*",
                        default=list(TARGET_ALGOS.keys()),
                        help="Targets to compile")
    args = parser.parse_args()

    for target in args.targets:
        print(f"\nCompiling {target}...")
        compile_target(target, args.output_dir)


if __name__ == "__main__":
    main()
