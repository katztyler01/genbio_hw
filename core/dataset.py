import pandas as pd
import pysam
import torch
from torch.utils.data import Dataset


class ATACCageSeqDataset(Dataset):
    def __init__(self, data_parquet: str, ref_genome: str):
        self.df = pd.read_parquet(data_parquet)
        self.ref_genome = pysam.FastaFile(ref_genome)
    
    def __len__(self):
        return self.df.shape[0]
    
    def __getitem__(self, idx):
        # NOTE: This returns the sequnce as a variable length string, not a tensor. Would tokenize here, chunk to context size in preprocessing
        row = self.df.iloc[idx]
        chrom = row["chrom"]
        start = int(row["start"])
        end = int(row["end"])
        assay = row["assay"]
        uniprot_ids = row["uniprot_ids"]
        
        seq = self._fetch_sequence(chrom, start, end)
        return seq, assay, uniprot_ids
    
    def _fetch_sequence(self, chrom: str, start: int, end: int):
        seq = self.ref_genome.fetch(chrom, start, end)
        return seq