import pandas as pd
import argparse
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge collected CSV and SMILES CSV")

    parser.add_argument(
        "--collected",
        type=str,
        required=True,
        help="Path to the collected CSV file"
    )
    parser.add_argument(
        "--smiles_csv",
        type=str,
        required=True,
        help="Path to the SMILES CSV file"
    )
    parser.add_argument(
        "--outdir",
        type=str,
        required=True,
        help="Output directory to save summary.csv"
    )
    parser.add_argument(
        "--smiles_col",
        type=str,
        default="SMILES",
        help="Column name in the SMILES CSV to merge on (default: 'smiles')"
    )

    args = parser.parse_args()

    df_sum = pd.read_csv(args.collected)
    df_smiles = pd.read_csv(args.smiles_csv)

    # Corrected: apply split to each row using .apply()
    df_sum['index'] = df_sum['Folder'].apply(lambda x: int(x.split('_')[-1]))

    # Merge on index
    df = pd.merge(df_sum, df_smiles[[args.smiles_col]], left_on='index', right_index=True, how='left')

    # Ensure output directory exists
    os.makedirs(args.outdir, exist_ok=True)
    outname = os.path.join(args.outdir, "summary.csv")

    df.to_csv(outname, index=False)
