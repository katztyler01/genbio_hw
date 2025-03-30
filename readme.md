# GenBio AI 
## Sequence Data Engineering Homework Assignment

### Tyler Katz

In `core/downloader.py`, an `ENCODEDownloader` class is defined to utilize real ENCODE API endpoints to search and download relevant files for ATAC-seq and CAGE experiments, as well as a human reference genome and annotation. 

Example Usage:
```
from core.downloader import ENCODEDownloader
downloader = ENCODEDownloader(
    output_dir="encode_data"
)
downloader.download_reference_genome()
downloader.find_experiments(
    assay="ATAC-seq",
    n_experiments=5
)
```

In `core/processor.py`, a `ENCODEProcessor` class is defined to process downloaded experiment files by extending peaks in either direction by 10 kb and merging them. Then, it finds genes annotated in these regions and utilizes the UniProt API to map Ensembl IDs to UniProt IDs for each isoform found. It saves the experiment as a Pandas DataFrame as a parquet file, where each row is an extended/merged peak region with additional experimental metadata.

Example Usage:
```
from core.processor import ENCODEProcessor

processor = ENCODEProcessor(
    output_dir="encode_data",
    ref_genome="encode_data/genomes/GRCh38.fa",
    db="encode_data/genomes/GRCh38_v24_annotation.sqlite",
    annotation="encode_data/genomes/GRCh38_v24.gtf",
    chrom_sizes="encode_data/genomes/GRCh38.chrom.sizes.tsv",
    ensembl_to_uniprot="encode_data/ensembl_to_uniprot.pkl",
    assembly"GRCh38"
)

processor.process_exp(("encode_data/experiments/ATAC-seq/ENCSR000RBT/ENCFF246XGY.bed", "encode_data/experiments/ATAC-seq/ENCSR000RBT/metadata.json" ))

```

