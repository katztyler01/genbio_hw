import pandas as pd
import os
import glob
from tqdm import tqdm
import argparse


def main(args):
    path = args.path
    parquet_files = glob.glob(os.path.join(path, "*.parquet"))

    dataframes = []
    for file in tqdm(parquet_files, desc="Loading parquet files", unit="file"):
        df = pd.read_parquet(file)
        dataframes.append(df)

    combined_df = pd.concat(dataframes, ignore_index=True)

    combined_df.to_parquet(args.out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Combine experiment parquet files into one."
    )
    parser.add_argument(
        "--path",
        type=str,
        default="encode_data/processed",
        help="Path to the directory containing parquet files.",
    )
    parser.add_argument(
        "--out",
        type=str,
        default="combined.parquet",
        help="Output file name for the combined parquet file.",
    )
    args = parser.parse_args()
    main(args)
