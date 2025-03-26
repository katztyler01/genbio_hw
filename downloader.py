import requests
import json
import os
from tqdm import tqdm
from typing import Optional, Dict

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

    def download_files(self, experiment: Dict, assay: str, file_types=None):
        exp_id = experiment["@id"]
        save_dir = f"{self.output_dir}/experiments/{assay}"
        save_dir = os.path.join(save_dir, experiment["accession"])
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        print(f"Experiment: {exp_id}")
        exp_url = f"https://www.encodeproject.org{exp_id}"
        response = requests.get(exp_url, headers={"accept": "application/json"})

        if response.status_code != 200:
            print(f"Error fetching experiment {exp_id}: {response.status_code}")
            return []
        exp_data = response.json()
        files = exp_data.get("files", [])

        if file_types:
            files = [f for f in files if f.get("file_type") in file_types]

        downloaded_files = []
        for file in files:
            file_url = f"https://www.encodeproject.org{file['href']}"
            file_name = os.path.join(
                save_dir, file["accession"] + "." + file["file_format"]
            )
            print(f"Downloading {file_url} to {file_name}")

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
            downloaded_files.append(file_name)
        return downloaded_files

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
            "limit": n_experiments if n_experiments else "all",
        }

        experiments = self.search_encode(params)
        for exp in experiments:
            self.download_files(
                exp, assay, file_types=["bed idr_peak", "bed idr_ranked_peak"]
            )
        return
