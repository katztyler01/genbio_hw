import argparse
import glob
import multiprocessing
import os
from functools import partial
from typing import Any, Dict, List, Optional, Tuple

from core.processor import ENCODEProcessor


def process_experiment(exp_args, args):
    bed_files, metadata_json = exp_args
    processor = ENCODEProcessor(
        output_dir=args.output_dir,
        ref_genome=args.ref_genome,
        annotation=args.annotation,
        db=args.db,
        chrom_sizes=args.chrom_sizes,
        ensembl_to_uniprot=args.ensembl_to_uniprot,
    )
    processor.process_exp((bed_files, metadata_json))


def main(args):
    assay_folder = os.path.join(f"{args.output_dir}/experiments", args.assay)
    experiment_data = []
    print(f"Finding ENCODE {args.assay} data...")
    for experiment_folder in os.listdir(assay_folder):
        experiment_path = os.path.join(assay_folder, experiment_folder)

        if not os.path.isdir(experiment_path):
            continue

        bed_files = glob.glob(os.path.join(experiment_path, "*.bed"))
        json_files = glob.glob(os.path.join(experiment_path, "*.json"))

        if bed_files and json_files:
            metadata_json = json_files[0]
            experiment_data.append((bed_files, metadata_json))

    n_workers = args.n_workers
    if n_workers is None:
        n_workers = multiprocessing.cpu_count()

    print(f"Processing {len(experiment_data)} experiments with {n_workers} workers...")
    process_fn = partial(process_experiment, args=args)
    with multiprocessing.Pool(processes=n_workers) as pool:
        pool.map(process_fn, experiment_data)

    print("Processing complete")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process ENCODE data.")
    parser.add_argument(
        "--output_dir", type=str, help="encode_data directory", default="encode_data"
    )
    parser.add_argument(
        "--ref_genome",
        type=str,
        help="Reference genome file",
        default="encode_data/genomes/GRCh38.fa",
    )
    parser.add_argument(
        "--annotation",
        type=str,
        help="Annotation file",
        default="encode_data/genomes/GRCh38_v24.gtf",
    )
    parser.add_argument(
        "--db",
        type=str,
        help="Annotation database file",
        default="encode_data/genomes/GRCh38_v24_annotation.sqlite",
    )
    parser.add_argument(
        "--chrom_sizes",
        type=str,
        help="Chromosome sizes file",
        default="encode_data/genomes/GRCh38.chrom.sizes.tsv",
    )
    parser.add_argument(
        "--ensembl_to_uniprot",
        type=str,
        help="Ensembl to UniProt mapping pkl file",
        default="encode_data/ensembl_to_uniprot.pkl",
    )
    parser.add_argument(
        "--assay", type=str, help="Assay type (ATAC-seq/CAGE)", default="ATAC-seq"
    )
    parser.add_argument("--n_workers", type=int, help="Number of workers", default=None)
    args = parser.parse_args()
    main(args)
