# TranscriptExplorer Pipeline

`TranscriptExplorer` is a comprehensive computational pipeline designed to process all human protein-coding genes, align their transcript isoforms, predict their membrane topologies (using **DeepTMHMM**) and subcellular localizations (using **DeepLoc2**), and compile their normalized expression data across cancer types (TCGA) and normal tissues (GTEx).

The generated database files are used to power the interactive plots in the [TranscriptExplorer Dash Application](https://github.com/apollis44/Transcripts_explorer_app).

---

## 📋 Table of Contents
1. [Project Structure](#-project-structure)
2. [Data Requirements & Pre-requisites](#-data-requirements--pre-requisites)
3. [Environment Setup](#-environment-setup)
   - [Method A: Standard Python venv & pip (Recommended)](#method-a-standard-python-venv--pip)
   - [Method B: Modern uv Package Manager (Fastest)](#method-b-modern-uv-package-manager)
4. [External Tool Configurations](#%EF%B8%8F-external-tool-configurations)
   - [DeepTMHMM Setup](#1-deeptmhmm-setup)
   - [MAFFT Alignment Setup](#2-mafft-alignment-setup)
5. [Running the Data Generation Pipeline](#-running-the-data-generation-pipeline)
   - [Step 1: Test the pipeline](#step-1-test-the-pipeline)
   - [Step 2: Run the full pipeline](#step-2-run-the-full-pipeline)
6. [Database Portability & Windows Conversion](#%EF%B8%8F-database-portability--windows-conversion)
7. [Deploying Data to the Dash App](#-deploying-data-to-the-dash-app)

---

## 📁 Project Structure

*   `generating_necessary_files.py`: The main production data generation script. It uses Python's `multiprocessing` (with `fork` context) to process all human proteins in parallel.
*   `generating_necessary_files_test.py`: A test version of the pipeline running on a small, hardcoded subset of 4 genes to verify system dependencies and correctness.
*   `convert_linux_db_to_win.py`: Converts Linux native shelve databases (GDBM format) to a universal format (`.bak`, `.dat` and `.dir` files) readable on Windows.
*   `scripts/`: Core pipeline logic helper modules:
    *   `scripts/topology.py`: Handles sequence retrieval, cache lookups, and membrane topology parsing.
    *   `scripts/expression.py`: Processes and merges TCGA/GTEx H5AD expression counts with phenotype metadata, generating boxplot statistics.
    *   `scripts/extract_sequences.py`: Contains API fetch methods for Ensembl (fetch transcripts, protein sequence, cDNA lengths).
    *   `scripts/files_genes/`: Contains input files, including the list of target genes, gene expression data and expression metadata files.
*   `deeploc-2.1.All.tar.gz`: Pre-packaged local archive for installing DeepLoc2.

---

## 📦 Data Requirements & Pre-requisites

Before running the pipeline, ensure the following inputs exist inside `scripts/files_genes/`:

1.  **Gene List (`genes.txt`)**: A CSV file containing Ensembl Gene IDs, names, and synonyms, downloaded from BioMart.
2.  **Explored Log (`already_explored_genes.pkl`)**: A Python pickle file containing a list of genes that have already been successfully processed. The pipeline uses this to skip completed genes and support resuming a run.
3.  **Expression Counts Matrix (`transcript_counts_sparse.h5ad`)**: A large H5AD file containing the sparse matrix of transcript expression counts across samples, downloaded from UCSC Xena browser using host: https://toil.xenahubs.net, cohort: TCGA TARGET GTEx and dataset: TcgaTargetGtex_rsem_isoform_tpm.
4.  **Phenotype Metadata (`TcgaTargetGTEX_phenotype.txt`)**: TSV file mapping sample IDs to their study cohort (TCGA vs GTEx) and detailed tissue/category type.

### Ensembl API Sequence Cache
To avoid rate limits and network latency from the Ensembl REST API, the pipeline runs offline/locally using a precomputed sequence cache.
*   You must have a persistent cache database located at: `.cache_ensembl/sequences_cache` (.dat, .dir, .bak) in the root of the repository.
*   If a gene Ensembl ID is not found in this cache, the pipeline will raise an exception.
*   **Building the cache**: If you need to build or update this cache from the Ensembl peptide FASTA file (`scripts/files_genes/Homo_sapiens.GRCh38.pep.all.fa`), run the utility script:
    ```bash
    uv run convert_fasta_to_cache.py
    ```
    This parses the FASTA file, groups identical protein sequences for transcripts under each gene to prevent duplicate runs, and builds the database.

---

## ⚙️ Environment Setup

> [!IMPORTANT]
> Python `3.9` is strictly required for compatibility with PyTorch 1.7.1, DeepTMHMM dependencies, and DeepLoc2.

If you have `uv` installed, setup is simplified and much faster using the provided `pyproject.toml`:
```bash
# Sync environment and build all dependencies
uv sync
# Activate the environment
source .venv/bin/activate
```

---

## 🛠️ External Tool Configurations

### 1. DeepTMHMM Setup
1.  Obtain the DeepTMHMM software license package from the [DTU Health Tech Website](https://services.healthtech.dtu.dk/services/DeepTMHMM-1.0/).
2.  Extract the downloaded `DeepTMHMM-1.0` archive into the root of this repository. Rename the folder to exactly `DeepTMHMM`.
3.  Copy the script `scripts/predict_api.py` into the `DeepTMHMM` directory:
    ```bash
    cp scripts/predict_api.py DeepTMHMM/
    ```
4.  Install the DeepTMHMM dependencies:
    ```bash
    uv pip install -r DeepTMHMM/requirements.txt
    ```

### 2. MAFFT Alignment Setup
The pipeline requires MAFFT for sequence alignment.
*   **Ubuntu/Debian:**
    ```bash
    sudo apt-get install -y mafft
    ```
*   **macOS (Homebrew):**
    ```bash
    brew install mafft
    ```
---

## 🚀 Running the Data Generation Pipeline

### Step 1: Test the pipeline
Before starting a full pipeline run, verify that PyTorch, DeepTMHMM models, DeepLoc2, and MAFFT are correctly configured by running the test script:
```bash
uv run generating_necessary_files_test.py
```
This script will test 4 genes (`ENSG00000177455`, `ENSG00000156738`, `ENSG00000003056`, `ENSG00000261857`)

### Step 2: Run the full pipeline
When the test script finishes successfully, run the production pipeline:
```bash
uv run generating_necessary_files.py
```
*   **Parallelism:** By default, it runs with 50 parallel worker processes using the `fork` context. Reduce this number if you are working on a small computer.
*   **Output:** Generates native shelve database files under `./files_for_plots/`.
*   **Resumability:** The pipeline automatically tracks progress. If it is stopped or interrupted, re-running the script will pick up where it left off using `already_explored_genes.pkl`. If you want to run from scratch, you should empty `already_explored_genes.pkl`.

---

## 🔄 Database Portability & Windows Conversion

Shelve databases generated on Linux (using the native `gdbm` module) are binary-incompatible with Windows. If you plan to host or run the Dash application on a Windows machine or on Render as it currently:

1.  Make sure the generated database files (`deeploc2_output`, `TCGA_GTEx_plotting_data`, `transcripts_to_isoforms_mapping`, `membrane_topology_objects`) are moved or copied to the folder specified in `convert_linux_db_to_win.py` (which defaults to `./files_for_plots`, or adjust the path in the script).
2.  Run the conversion script:
    ```bash
    uv run convert_linux_db_to_win.py
    ```

---

## 🚀 Deploying Data to the Dash App

To deploy the generated databases to the interactive web application, copy them to the Dash application's data directory:
```bash
# For Linux-native deployment
cp -r files_for_plots/* /path/to/Transcript_explorer_app/files_for_plots/

# For Windows-compatible deployment (after running convert_linux_db_to_win.py)
cp app/files_for_plots/*.dat /path/to/Transcript_explorer_app/files_for_plots/
cp app/files_for_plots/*.dir /path/to/Transcript_explorer_app/files_for_plots/
cp app/files_for_plots/*.bak /path/to/Transcript_explorer_app/files_for_plots/
```
The Dash application will read these files to serve plots interactively!