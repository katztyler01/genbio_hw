import copy
import json
import multiprocessing
import os
import pickle
from collections import defaultdict
from functools import partial
from typing import Any, Dict, List, Optional, Tuple

import gffutils
import pandas as pd
from pybedtools import BedTool
from tqdm import tqdm

from uniprot import UniprotAPI


class ENCODEProcessor:
    def __init__(
        self,
        output_dir: str,
        ref_genome: str,
        db: str,
        annotation: str,
        chrom_sizes: str,
        ensembl_to_uniprot: Optional[str] = None,
        assembly: Optional[str] = "GRCh38",
    ) -> None:
        """
        Initializes the ENCODEProcessor with the given parameters.
        Args:
            output_dir (str): The directory to save the output files.
            ref_genome (str): The reference genome file.
            db (str): The database file for gffutils.
            annotation (str): The annotation file (GTF).
            chrom_sizes (str): The chromosome sizes file.
            ensembl_to_uniprot (str, optional): Path to the Ensembl to UniProt mapping pkl file. Defaults to None.
            assembly (str, optional): Genome assembly version. Defaults to "GRCh38".
        """
        self.output_dir = output_dir
        self.annotation = annotation
        self.db = db
        self.annotation = annotation
        self.ref_genome = ref_genome
        self.chrom_sizes = chrom_sizes
        self.assembly = assembly

        self.gene_db = self.get_gene_db(gtf_file=annotation, db_file=db)

        self.genes_bed = BedTool(
            [
                (gene.seqid, gene.start - 1, gene.end, gene.id, 0, gene.strand)
                for gene in self.gene_db.features_of_type("gene")
            ]
        )

        self.transcript_to_gene, self.gene_to_transcripts = (
            self.get_transcript_gene_mapping()
        )
        self.ensembl_to_uniprot_pkl = ensembl_to_uniprot

    def get_gene_db(
        self, gtf_file="encode_data/genomes/GRCh38_v24.gtf.gz", db_file=None
    ) -> gffutils.FeatureDB:
        """
        Creates or loads a gffutils database from the provided GTF file.
        Args:
            gtf_file (str): Path to the GTF file.
            db_file (str, optional): Path to the database file. Defaults to None.
        Returns:
            gffutils.FeatureDB: The gffutils database object.
        """
        if db_file is None:
            db_file = ":memory:"
        elif os.path.exists(db_file):
            try:
                print(f"Loading existing database from {db_file}...")
                db = gffutils.FeatureDB(db_file)
                print("Database loaded successfully")
                return db
            except Exception as e:
                print(f"Error loading database: {e}")

        print("Creating new database...")
        db = gffutils.create_db(
            gtf_file,
            db_file,
            force=True,
            keep_order=True,
            merge_strategy="merge",
            sort_attribute_values=True,
            disable_infer_genes=True,
            disable_infer_transcripts=True,
        )

        return db

    def get_transcript_gene_mapping(
        self,
    ) -> Tuple[Dict[str, str], Dict[str, List[str]]]:
        """
        Creates a mapping of transcript IDs to gene IDs and vice versa.
        Returns:
            (Tuple[Dict[str, str], Dict[str, List[str]]]):
                - transcript_to_gene: Mapping of transcript IDs to gene IDs.
                - gene_to_transcripts: Mapping of gene IDs to a list of transcript IDs.
        """
        transcript_to_gene = {}
        gene_to_transcripts = defaultdict(list)
        for transcript in self.gene_db.features_of_type("transcript"):
            gene_type = transcript["gene_type"][0]
            transcript_type = transcript["transcript_type"][0]
            if gene_type != "protein_coding" or transcript_type != "protein_coding":
                continue
            transcript_id = transcript["transcript_id"][0]
            gene_id = transcript["gene_id"][0]
            transcript_to_gene[transcript_id] = gene_id
            gene_to_transcripts[gene_id].append(transcript_id)

        return transcript_to_gene, gene_to_transcripts

    def map_ensembl_to_uniprot(
        self, transcript_ids: List[str], gene_ids: Optional[List[str]] = None
    ) -> Dict[str, str]:
        """
        Maps Ensembl IDs to UniProt IDs using the UniProt API.
        Args:
            transcript_ids (List[str]): List of Ensembl transcript IDs.
            gene_ids (List[str], optional): List of Ensembl gene IDs. Defaults to None.
        Returns:
            ensembl_to_uniprot (Dict[str, str]): Mapping of Ensembl IDs to UniProt IDs.
        """
        # NOTE: code adapted from UniProt API docs https://www.uniprot.org/help/id_mapping
        uniprot_api = UniprotAPI()
        ensembl_to_uniprot = {}

        if self.ensembl_to_uniprot_pkl and os.path.exists(self.ensembl_to_uniprot_pkl):
            with open(self.ensembl_to_uniprot_pkl, "rb") as f:
                ensembl_to_uniprot = pickle.load(f)

            transcript_ids = [
                tid
                for tid in transcript_ids
                if tid.split(".")[0] not in ensembl_to_uniprot
            ]
            n_gene_ids = 0
            if gene_ids:
                gene_ids = [
                    gid
                    for gid in gene_ids
                    if gid.split(".")[0] not in ensembl_to_uniprot
                ]
                n_gene_ids = len(gene_ids)

            total_ids = len(transcript_ids) + n_gene_ids
            if total_ids == 0:
                return ensembl_to_uniprot

        batch_size = 50000
        for i in range(0, len(transcript_ids), batch_size):
            batch = transcript_ids[i : i + batch_size]

            base_ids = []
            for tid in batch:
                base_ids.append(tid.split(".")[0])

            job_id = uniprot_api.submit_id_mapping(
                from_db="Ensembl_Transcript", to_db="UniProtKB", ids=base_ids
            )
            if uniprot_api.check_id_mapping_results_ready(job_id):
                link = uniprot_api.get_id_mapping_results_link(job_id)
                results = uniprot_api.get_id_mapping_results_search(link)

            for result in results["results"]:
                transcript_id = result["from"]
                uniprot_id = result["to"]["primaryAccession"]
                ensembl_to_uniprot[transcript_id] = uniprot_id

        if gene_ids:
            for i in range(0, len(gene_ids), batch_size):
                batch = gene_ids[i : i + batch_size]

                base_ids = []
                for gid in batch:
                    base_ids.append(gid.split(".")[0])

                job_id = uniprot_api.submit_id_mapping(
                    from_db="Ensembl", to_db="UniProtKB", ids=base_ids
                )
                if uniprot_api.check_id_mapping_results_ready(job_id):
                    link = uniprot_api.get_id_mapping_results_link(job_id)
                    results = uniprot_api.get_id_mapping_results_search(link)

                for result in results["results"]:
                    gene_id = result["from"]
                    uniprot_id = result["to"]["primaryAccession"]
                    ensembl_to_uniprot[gene_id] = uniprot_id

        if self.ensembl_to_uniprot_pkl and not os.path.exists(
            self.ensembl_to_uniprot_pkl
        ):
            with open(self.ensembl_to_uniprot_pkl, "wb") as f:
                pickle.dump(ensembl_to_uniprot, f)

        return ensembl_to_uniprot

    def process_all_exps(self, bed_files: List[str], output: str, assay: str) -> None:
        """
        Processes all experiments for a given assay together
        and saves a combined BED file and Pandas DataFrame as a parquet file.
        Args:
            bed_files (List[str]): List of paths to BED files.
            output (str): Path to save the output parquet file.
            assay (str): The type of assay (ATAC-seq, CAGE).
        """
        peaks = BedTool(bed_files[0])
        for i, bed_file in tqdm(
            enumerate(bed_files),
            total=len(bed_files),
            desc="Processing bed files",
            unit="file",
        ):
            if i == 0:
                continue
            bed = BedTool(bed_file)
            peaks = peaks.cat(bed)

        peaks = peaks.slop(b=10000, g=self.chrom_sizes)
        peaks = peaks.sort().merge()
        peaks.saveas(f"{self.output_dir}/{assay}_all.bed")
        peak_gene_intersect = peaks.intersect(self.genes_bed, loj=True)

        peak_to_genes = defaultdict(set)
        all_gene_ids = set()

        for peak in peak_gene_intersect:
            chrom = peak[0]
            peak_start = int(peak[1])
            peak_end = int(peak[2])
            gene_id = peak[6]
            if gene_id == ".":
                continue
            peak_to_genes[(chrom, peak_start, peak_end)].add(gene_id)
            all_gene_ids.add(gene_id)

        transcript_ids = set()
        for gene_id in all_gene_ids:
            if gene_id in self.gene_to_transcripts:
                transcript_ids.update(self.gene_to_transcripts[gene_id])

        # ensembl_to_uniprot = self.map_ensembl_to_uniprot(list(transcript_ids), list(all_gene_ids))
        ensembl_to_uniprot = self.map_ensembl_to_uniprot(list(transcript_ids))

        samples = []
        for peak_id, gene_ids in peak_to_genes.items():
            sample = {}
            sample["chrom"] = peak_id[0]
            sample["start"] = peak_id[1]
            sample["end"] = peak_id[2]
            sample["genes"] = list(gene_ids)
            sample["assay"] = assay
            sample["assembly"] = self.assembly

            uniprot_ids = set()
            for gene_id in gene_ids:
                stable_gene_id = gene_id.split(".")[0]
                uniprot_id = ensembl_to_uniprot.get(stable_gene_id)
                if uniprot_id:
                    uniprot_ids.add(uniprot_id)
                if gene_id in self.gene_to_transcripts:
                    transcript_ids = self.gene_to_transcripts[gene_id]
                    for transcript_id in transcript_ids:
                        stable_transcript_id = transcript_id.split(".")[0]
                        if stable_transcript_id in ensembl_to_uniprot:
                            uniprot_id = ensembl_to_uniprot[stable_transcript_id]
                            uniprot_ids.add(uniprot_id)

            sample["uniprot_ids"] = list(uniprot_ids)
            samples.append(sample)

        df = pd.DataFrame(samples)
        df.to_parquet(output)
        return

    def process_exp(self, args: Tuple[List[str], str]) -> None:
        """
        Processes a single experiment by extending peaks and merging BED files and saves processed experiment as Pandas DataFrame.
        Args:
            args (Tuple[List[str], str]): Tuple containing a list of BED files and the metadata JSON file.
        """
        bed_files, metadata_file = args
        with open(metadata_file, "r") as file:
            metadata = json.load(file)

        assay = metadata.get("assay_title")
        if assay is None:
            print(f"Assay title not found in metadata file {metadata_file}")
            return

        df_output = os.path.join(
            self.output_dir, f"processed_{assay}/{metadata['accession']}.parquet"
        )

        if os.path.exists(df_output):
            print(f"Already processed experiment {metadata['accession']}")
            return

        if not os.path.exists(os.path.dirname(df_output)):
            os.makedirs(os.path.dirname(df_output))

        annotations = {}
        annotations["cell_slims"] = metadata["biosample_ontology"].get("cell_slims", [])
        annotations["organ_slims"] = metadata["biosample_ontology"].get(
            "organ_slims", []
        )
        annotations["system_slims"] = metadata["biosample_ontology"].get(
            "system_slims", []
        )
        annotations["biosample_class"] = metadata["biosample_ontology"].get(
            "classification"
        )
        annotations["assembly"] = self.assembly
        annotations["assay"] = assay

        peaks = BedTool(bed_files[0])
        for i, bed_file in enumerate(bed_files):
            if i == 0:
                continue
            bed = BedTool(bed_file)
            peaks = peaks.cat(bed, postmerge=False)

        peaks = peaks.slop(b=10000, g=self.chrom_sizes)
        peaks = peaks.sort().merge()
        peak_gene_intersect = peaks.intersect(self.genes_bed, loj=True)

        peak_to_genes = defaultdict(set)
        all_gene_ids = set()

        for peak in peak_gene_intersect:
            chrom = peak[0]
            peak_start = int(peak[1])
            peak_end = int(peak[2])
            gene_id = peak[6]
            if gene_id == ".":
                continue
            peak_to_genes[(chrom, peak_start, peak_end)].add(gene_id)
            all_gene_ids.add(gene_id)

        transcript_ids = set()
        for gene_id in all_gene_ids:
            if gene_id in self.gene_to_transcripts:
                transcript_ids.update(self.gene_to_transcripts[gene_id])

        # ensembl_to_uniprot = self.map_ensembl_to_uniprot(list(transcript_ids), list(all_gene_ids))
        ensembl_to_uniprot = self.map_ensembl_to_uniprot(list(transcript_ids))

        samples = []
        for peak_id, gene_ids in peak_to_genes.items():
            sample = copy.deepcopy(annotations)
            sample["chrom"] = peak_id[0]
            sample["start"] = peak_id[1]
            sample["end"] = peak_id[2]
            sample["genes"] = list(gene_ids)

            uniprot_ids = set()
            for gene_id in gene_ids:
                stable_gene_id = gene_id.split(".")[0]
                uniprot_id = ensembl_to_uniprot.get(stable_gene_id)
                if uniprot_id:
                    uniprot_ids.add(uniprot_id)
                if gene_id in self.gene_to_transcripts:
                    transcript_ids = self.gene_to_transcripts[gene_id]
                    for transcript_id in transcript_ids:
                        stable_transcript_id = transcript_id.split(".")[0]
                        if stable_transcript_id in ensembl_to_uniprot:
                            uniprot_id = ensembl_to_uniprot[stable_transcript_id]
                            uniprot_ids.add(uniprot_id)

            sample["uniprot_ids"] = list(uniprot_ids)
            samples.append(sample)

        df = pd.DataFrame(samples)
        df.to_parquet(df_output)

        # return samples
        return
