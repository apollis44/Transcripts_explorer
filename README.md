# TranscriptExplorer

TranscriptExplorer is a comprehensive pipeline and visualization tool designed to analyze and explore the membrane topology and expression of different protein isoforms across various cancer types (using TCGA data) and normal tissues (using GTEx data) for all human proteins.

## Project Structure

- `generating_necessary_files.py`: The main data generation script. It processes a list of genes, retrieves transcript data, aligns sequences, predicts membrane topologies (using DeepTMHMM via Biolib), and compiles expression data. This script only has to be ran once on all human proteins.
- `scripts/`: Contains Python modules for extracting sequences, predicting topology, and processing expression data.
- `app/`: Contains a Plotly Dash web application for visualizing the generated data interactively.
- `files/`: Output directory where intermediate sequences, functional predictions, and alignments are saved temporarily.

## Setup

1. Clone the repository.
2. Download DeepTMHMM from [here](https://services.healthtech.dtu.dk/services/DeepTMHMM-1.0/) and put the unzipped folder in the Transcript_explorer folder.
3. Install the required Python dependencies:
   ```bash
   cd Transcript_explorer
   sudo apt-get install python3.9 python3.9-dev python3.9-venv libhdf5-dev
   python3.9 -m venv ../TranscriptExplorer_env
   source ../TranscriptExplorer_env/bin/activate
   pip install wheel Cython==0.29.37 pkgconfig==1.5.5
   pip install https://download.pytorch.org/whl/cpu/torch-1.7.1%2Bcpu-cp39-cp39-linux_x86_64.whl#sha256=445c5ff49964a3cdd8170a20a7371e3691412e2e4a7005f9c89c485ab47e8609
   pip install -r DeepTMHMM/requirements.txt
   pip install -r requirements.txt
   pip install "pytorch_lightning<1.5.0" "torchmetrics<0.8.0" "transformers<4.20.0"
   pip install deeploc-2.1.All.tar.gz --upgrade-strategy only-if-needed --no-cache-dir
   ```
4. Copy-paste the scripts/predict_api.py file into the DeepTMHMM folder.
5. Install MAFFT (Multiple-Sequence Alignment program):
   - On Ubuntu/Debian: `sudo apt-get install mafft`
   - On macOS (Homebrew): `brew install mafft`
   - Or download from the [MAFFT website](https://mafft.cbrc.jp/alignment/software/).

6. Test the pipeline:
   ```bash
   python3 generating_necessary_files_test.py
   ```

7. Generate the data mapping and objects:
   ```bash
   python3 generating_necessary_files.py
   ```
   *Note: Ensure you configure your Biolib API token and use your own email in the script.*

8. To view the results interactively, navigate to the `app/` directory and run the Dash server. (See `app/README.md` for more details).