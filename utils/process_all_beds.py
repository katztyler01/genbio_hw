from core.processor import ENCODEProcessor
import os
import argparse
import glob


def main(args):
    processor = ENCODEProcessor(
        output_dir=args.output_dir,
        ref_genome=args.ref_genome,
        annotation=args.annotation,
        db=args.db,
        chrom_sizes=args.chrom_sizes,
        ensembl_to_uniprot=args.ensembl_to_uniprot,
    )
    all_bed_files = []
    for experiment_folder in os.listdir(args.assay_folder):
        experiment_path = os.path.join(args.assay_folder, experiment_folder)

        if not os.path.isdir(experiment_path):
            continue

        bed_files = glob.glob(os.path.join(experiment_path, "*.bed"))

        if bed_files:
            all_bed_files.extend(bed_files)

    processor.process_all_exps(
        bed_files=all_bed_files,
        output=f"{args.output_dir}/{args.assay}_all.parquet",
        assay=args.assay,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process all experiments for a given assay together."
    )
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
    parser.add_argument(
        "--assay_folder",
        type=str,
        help="Folder containing the assay data",
        default="encode_data/experiments/ATAC-seq",
    )
    args = parser.parse_args()
    main(args)
