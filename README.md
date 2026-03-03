# TranscriptExplorer

TranscriptExplorer is a comprehensive pipeline and visualization tool designed to analyze and explore the membrane topology and expression of different protein isoforms across various cancer types (using TCGA data) and normal tissues (using GTEx data) for all human proteins.

## Project Structure

- `generating_necessary_files.py`: The main data generation script. It processes a list of genes, retrieves transcript data, aligns sequences, predicts membrane topologies (using DeepTMHMM via Biolib), and compiles expression data. This script only has to be ran once on all human proteins.
- `scripts/`: Contains Python modules for extracting sequences, predicting topology, and processing expression data.
- `app/`: Contains a Plotly Dash web application for visualizing the generated data interactively.
- `files/`: Output directory where intermediate sequences, functional predictions, and alignments are saved.

## Setup

1. Clone the repository.
2. Install the required Python dependencies:
   ```bash
   pip install -r requirements.txt
   pip install deeploc-2.1.All.tar.gz
   ```
3. Install MAFFT (Multiple-Sequence Alignment program):
   - **Windows:** MAFFT is already included in this repository and accessible.
   - **Linux/macOS:** You must install MAFFT on your system.
     - On Ubuntu/Debian: `sudo apt-get install mafft`
     - On macOS (Homebrew): `brew install mafft`
     - Or download from the [MAFFT website](https://mafft.cbrc.jp/alignment/software/).

4. Generate the data mapping and objects:
   ```bash
   python generating_necessary_files.py
   ```
   *Note: Ensure you configure your Biolib API token and use your own email in the script.*

5. To view the results interactively, navigate to the `app/` directory and run the Dash server. (See `app/README.md` for more details).