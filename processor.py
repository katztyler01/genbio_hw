import gffutils
import os
from collections import defaultdict
from downloader import ENCODEDownloader
from pybedtools import BedTool
from typing import Dict, List, Tuple

class ENCODEProcessor:
    def __init__(self, output_dir: str, db_file: str) -> None:
        self.output_dir = output_dir
        self.downloader = ENCODEDownloader(output_dir)
        self.gene_db = self.downloader.get_gene_db(db_file=db_file)
    
    def get_transcript_gene_mapping(self) -> Tuple[Dict[str, str], Dict[str, List[str]]]:
        transcript_to_gene = {}
        gene_to_transcripts = defaultdict(list)
        for transcript in self.gene_db.features_of_type('transcript'):
            gene_type = transcript['gene_type'][0]
            if gene_type != 'protein_coding':
                continue
            transcript_id = transcript['transcript_id'][0]
            gene_id = transcript['gene_id'][0]
            transcript_to_gene[transcript_id] = gene_id
            gene_to_transcripts[gene_id].append(transcript_id)
            
        self.transcript_to_gene = transcript_to_gene
        self.gene_to_transcripts = gene_to_transcripts
        
        return transcript_to_gene, gene_to_transcripts
    
    def process_bed(self, bed_file: str, genome: str) -> None:
        peaks = BedTool(bed_file)
        peaks = peaks.slop(b=10000, genome=genome)
        peaks = peaks.sort().merge()
