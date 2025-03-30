import argparse

import pandas as pd
import pysam


def load_reference_genome(fasta_file):
    ref = pysam.FastaFile(fasta_file)
    return ref


def calculate_gc_content(sequence):
    if not sequence:
        return None

    g_count = sequence.upper().count("G")
    c_count = sequence.upper().count("C")

    return (g_count + c_count) / len(sequence) * 100


def add_gc_content(df, ref):
    def get_gc_content(row):
        chrom = row["chrom"]
        start = int(row["start"])
        end = int(row["end"])
        seq = ref.fetch(chrom, start, end)
        return calculate_gc_content(seq)

    df["gc_content"] = df.apply(get_gc_content, axis=1)
    return df


def main(args):
    print(f"Loading parquet file: {args.parquet}")
    df = pd.read_parquet(args.parquet)
    print(f"Loading reference genome: {args.fasta}")
    ref_genome = load_reference_genome(args.fasta)

    print(f"Calculating GC content for {len(df)} regions")
    df_gc = add_gc_content(df, ref_genome)
    df_gc.to_parquet(args.output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate GC content from a parquet file."
    )
    parser.add_argument("--parquet", help="Input parquet file")
    parser.add_argument("--fasta", help="Reference genome in FASTA format")
    parser.add_argument("--output", help="Output parquet file with GC content")

    args = parser.parse_args()

    main(args)
