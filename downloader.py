import requests
import json
import os
from tqdm import tqdm
from typing import Optional, Dict
from pybedtools import BedTool
import gffutils
import gzip
import shutil

ENCODE_BASE_URL = "https://www.encodeproject.org"
ENCODE_SEARCH_URL = "https://www.encodeproject.org/search/"
ENCODE_FILE_URL = "https://www.encodeproject.org/files/"


class ENCODEDownloader:
    def __init__(self, output_dir: str = "./encode_data"):
        self.output_dir = output_dir

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

    def search_encode(self, query_params):
        headers = {"accept": "application/json"}
        response = requests.get(ENCODE_SEARCH_URL, params=query_params, headers=headers)

        if response.status_code != 200:
            print(f"Error searching ENCODE: {response.status_code}")
            response.raise_for_status()

        data = response.json()
        print(f"Found {len(data['@graph'])} entries matching the query")

        return data["@graph"]
    
    def download_file(self, file_url, file_name):
        with requests.get(file_url, stream=True) as r:
            r.raise_for_status()
            with open(file_name, "wb") as f, tqdm(
                desc=file_name,
                total=int(r.headers.get("content-length", 0)),
                unit="iB",
                unit_scale=True,
                unit_divisor=1024,
            ) as progress_bar:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
                    progress_bar.update(len(chunk))
    

    def download_files(self, experiment: Dict, assay: str):
        exp_id = experiment["@id"]
        save_dir = f"{self.output_dir}/experiments/{assay}"
        save_dir = os.path.join(save_dir, experiment["accession"])
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        print(f"Experiment: {exp_id}")
        exp_url = f"{ENCODE_BASE_URL}{exp_id}"
        response = requests.get(exp_url, headers={"accept": "application/json"})

        if response.status_code != 200:
            print(f"Error fetching experiment {exp_id}: {response.status_code}")
            return []
        exp_data = response.json()
        files = exp_data.get("files", [])

        
        idr_thresholded = []
        idr_ranked = []
        for f in files:
            if f.get("output_type") in ["IDR thresholded peaks"] and f.get("file_type").startswith("bed"):
                idr_thresholded.append(f)
            elif f.get("output_type") in ["IDR ranked peaks"] and f.get("file_type").startswith("bed"):
                idr_ranked.append(f)
                
        if len(idr_thresholded) > 0:
            files = idr_thresholded
        elif len(idr_ranked) < 1:
            print(f"No IDR thresholded or ranked peaks found for {exp_id}")
        else:
            files = idr_ranked

        downloaded_files = []
        for file in files:
            file_url = f"{ENCODE_BASE_URL}{file['href']}"
            file_name = os.path.join(
                save_dir, file["accession"] + "." + file["file_format"]
            )
            print(f"Downloading {file_url} to {file_name}")
            self.download_file(file_url, file_name)
            downloaded_files.append(file_name)
            
        return downloaded_files
    
    def download_reference_genome(self):
        save_dir = f"{self.output_dir}/genomes"
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
            
        grch38v24gtf_url = "https://www.encodeproject.org/files/ENCFF343XGU/@@download/ENCFF343XGU.gtf.gz"
        grch38fa_url = "https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz"
        grch38_chrom_sizes_url = "https://www.encodeproject.org/files/GRCh38_EBV.chrom.sizes/@@download/GRCh38_EBV.chrom.sizes.tsv"
        
        gtf_file = f"{save_dir}/GRCh38_v24.gtf.gz"
        fa_file = f"{save_dir}/GRCh38.fa.gz"
        chrom_sizes_file = f"{save_dir}/GRCh38.chrom.sizes.tsv"
        
        if not os.path.exists(gtf_file):
            self.download_file(file_url=grch38v24gtf_url, file_name=gtf_file)
        
        if not os.path.exists(fa_file):
            self.download_file(file_url=grch38fa_url, file_name=fa_file)
        
        if not os.path.exists(chrom_sizes_file):
            self.download_file(file_url=grch38_chrom_sizes_url, file_name=chrom_sizes_file)
        
        if not os.path.exists(gtf_file[:-3]):
            with gzip.open(gtf_file, 'rb') as f_in:
                    with open(gtf_file[:-3], 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                        
        if not os.path.exists(fa_file[:-3]):
            with gzip.open(fa_file, 'rb') as f_in:
                    with open(fa_file[:-3], 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
    
    def get_gene_db(self, gtf_file="encode_data/genomes/GRCh38_v24.gtf.gz", db_file=None):
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
            merge_strategy='merge',
            sort_attribute_values=True,
            disable_infer_genes=True,
            disable_infer_transcripts=True
        )
        
        return db
        
    def download_experiments(self, assay: str, n_experiments: Optional[int] = None):
        assert assay in ["ATAC-seq", "CAGE"], (
            "Assay must be either 'CAGE' or 'ATAC-seq'"
        )

        params = {
            "type": "Experiment",
            "assay_title": assay,
            "replicates.library.biosample.donor.organism.scientific_name": "Homo sapiens",
            "frame": "object",
            "format": "json",
            "audit.ERROR.category!": "extremely low read depth",
            "limit": n_experiments if n_experiments else "all",
        }
        if assay == "ATAC-seq":
            params["audit.NOT_COMPLIANT.category!"] = "low FRiP score"

        experiments = self.search_encode(params)
        for exp in experiments:
            self.download_files(
                exp, assay
            )
        return