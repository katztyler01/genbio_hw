import gzip
import json
import os
import shutil
from typing import Any, Dict, List, Optional

import gffutils
import requests
from pybedtools import BedTool
from tqdm import tqdm

ENCODE_BASE_URL = "https://www.encodeproject.org"
ENCODE_SEARCH_URL = "https://www.encodeproject.org/search/"
ENCODE_FILE_URL = "https://www.encodeproject.org/files/"


class ENCODEDownloader:
    def __init__(self, output_dir: str = "./encode_data", assembly="GRCh38"):
        self.output_dir = output_dir

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        self.assembly = assembly

    def search_encode(self, query_params: Dict[str, Any]) -> Dict:
        """
        Uses ENCODE API to search using the provided query parameters.

        Args:
            query_params (Dict[str, Any]): The query parameters for the search.

        Returns:
            Dict: The search results in JSON format.
        """
        headers = {"accept": "application/json"}
        response = requests.get(ENCODE_SEARCH_URL, params=query_params, headers=headers)

        if response.status_code != 200:
            print(f"Error searching ENCODE: {response.status_code}")
            response.raise_for_status()

        data = response.json()
        print(f"Found {len(data['@graph'])} entries matching the query")

        return data["@graph"]

    def download_file(self, file_url: str, file_name: str) -> None:
        """
        Downloads a file from ENCODE using API.

        Args:
            file_url (str): The URL of the file to download.
            file_name (str): The name to save the downloaded file as.
        """
        with requests.get(file_url, stream=True) as r:
            r.raise_for_status()
            with (
                open(file_name, "wb") as f,
                tqdm(
                    desc=file_name,
                    total=int(r.headers.get("content-length", 0)),
                    unit="iB",
                    unit_scale=True,
                    unit_divisor=1024,
                ) as progress_bar,
            ):
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
                    progress_bar.update(len(chunk))

    def download_experiment(self, experiment: Dict, assay: str) -> Optional[List[str]]:
        """
        Retreives info about a specific experiment and downloads the relevant files.
        Args:
            experiment (Dict): The experiment data from ENCODE API as JSON.
            assay (str): The type of assay (ATAC-seq, CAGE).

        """
        exp_id = experiment["@id"]
        save_dir = f"{self.output_dir}/experiments/{assay}"
        save_dir = os.path.join(save_dir, experiment["accession"])

        print(f"Experiment: {exp_id}")
        exp_url = f"{ENCODE_BASE_URL}{exp_id}"
        response = requests.get(exp_url, headers={"accept": "application/json"})

        if response.status_code != 200:
            print(f"Error fetching experiment {exp_id}: {response.status_code}")
            return []
        exp_data = response.json()
        files = exp_data.get("files", [])

        if assay == "ATAC-seq":
            idr_thresholded = []
            idr_ranked = []
            for f in files:
                if f.get("output_type") in ["IDR thresholded peaks"] and f.get(
                    "file_type"
                ).startswith("bed"):
                    idr_thresholded.append(f)
                elif f.get("output_type") in ["IDR ranked peaks"] and f.get(
                    "file_type"
                ).startswith("bed"):
                    idr_ranked.append(f)

            if len(idr_thresholded) > 0:
                files = idr_thresholded
            elif len(idr_ranked) < 1:
                print(f"No IDR thresholded or ranked peaks found for {exp_id}")
                return None
            else:
                files = idr_ranked
        else:
            idr_peak = []
            for f in files:
                assembly = f.get("assembly")
                if assembly and assembly != self.assembly:
                    continue
                if f.get("file_type") in ["bed idr_peak"]:
                    idr_peak.append(f)
            if len(idr_peak) > 0:
                files = idr_peak
            else:
                print(f"No IDR peak files found for {exp_id}")
                return None

        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        metadata_file = os.path.join(save_dir, "metadata.json")
        with open(metadata_file, "w") as file:
            json.dump(experiment, file, indent=4)

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

    def download_reference_genome(self) -> None:
        """
        Downloads the reference genome and annotation files for GRCh38 from ENCODE.
        """
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
            self.download_file(
                file_url=grch38_chrom_sizes_url, file_name=chrom_sizes_file
            )

        if not os.path.exists(gtf_file[:-3]):
            with gzip.open(gtf_file, "rb") as f_in:
                with open(gtf_file[:-3], "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)

        if not os.path.exists(fa_file[:-3]):
            with gzip.open(fa_file, "rb") as f_in:
                with open(fa_file[:-3], "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)

    def find_experiments(self, assay: str, n_experiments: Optional[int] = None) -> None:
        """
        Finds and downloads experiments from ENCODE based on the specified assay type.

        Args:
            assay (str): The type of assay (ATAC-seq, CAGE).
            n_experiments (Optional[int]): The number of experiments to download. If None, all experiments will be downloaded.
        """
        assert assay in ["ATAC-seq", "CAGE"], (
            "Assay must be either 'CAGE' or 'ATAC-seq'"
        )

        params = {
            "type": "Experiment",
            "assay_title": assay,
            "replicates.library.biosample.donor.organism.scientific_name": "Homo sapiens",
            "frame": "embedded",
            "format": "json",
            "assembly": "GRCh38",
            "audit.ERROR.category!": "extremely low read depth",
            "limit": n_experiments if n_experiments else "all",
        }

        if assay == "ATAC-seq":
            params["audit.NOT_COMPLIANT.category!"] = "low FRiP score"

        experiments = self.search_encode(params)
        for exp in experiments:
            self.download_experiment(exp, assay)
        return
