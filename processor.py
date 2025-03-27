import gffutils
import os
import requests
from collections import defaultdict
from downloader import ENCODEDownloader
from pybedtools import BedTool
from xml.etree import ElementTree
import json
import time
import re
from urllib.parse import urlparse, parse_qs, urlencode
import zlib
from requests.adapters import HTTPAdapter, Retry
from typing import Dict, List, Tuple
from uniprot import UniprotAPI


class ENCODEProcessor:
    def __init__(self, output_dir: str, db_file: str) -> None:
        self.output_dir = output_dir
        self.downloader = ENCODEDownloader(output_dir)
        self.gene_db = self.downloader.get_gene_db(db_file=db_file)

    def get_transcript_gene_mapping(
        self,
    ) -> Tuple[Dict[str, str], Dict[str, List[str]]]:
        transcript_to_gene = {}
        gene_to_transcripts = defaultdict(list)
        for transcript in self.gene_db.features_of_type("transcript"):
            gene_type = transcript["gene_type"][0]
            if gene_type != "protein_coding":
                continue
            transcript_id = transcript["transcript_id"][0]
            gene_id = transcript["gene_id"][0]
            transcript_to_gene[transcript_id] = gene_id
            gene_to_transcripts[gene_id].append(transcript_id)

        self.transcript_to_gene = transcript_to_gene
        self.gene_to_transcripts = gene_to_transcripts

        return transcript_to_gene, gene_to_transcripts

    def map_transcripts_to_uniprot(self, transcript_ids):
        # NOTE: code adapted from UniProt API docs https://www.uniprot.org/help/id_mapping
        uniprot_api = UniprotAPI()
        transcript_to_uniprot = {}

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
                transcript_to_uniprot[transcript_id] = uniprot_id

        return transcript_to_uniprot

    def process_bed(self, bed_file: str, genome: str) -> None:
        peaks = BedTool(bed_file)
        peaks = peaks.slop(b=10000, genome=genome)
        peaks = peaks.sort().merge()
