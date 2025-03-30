from core.downloader import ENCODEDownloader
import argparse


def main(args):
    downloader = ENCODEDownloader(output_dir=args.output_dir)

    downloader.find_experiments(assay=args.assay, n_experiments=args.n_experiments)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download ENCODE data")
    parser.add_argument(
        "--output_dir",
        type=str,
        default="./encode_data",
        help="Directory to save downloaded data",
    )
    parser.add_argument(
        "--assay",
        type=str,
        default="ATAC-seq",
        help="Assay type to download (ATAC-seq/CAGE)",
    )
    parser.add_argument(
        "--n_experiments",
        type=int,
        default=None,
        help="Numer of experiments to download",
    )
    args = parser.parse_args()
    main(args)
