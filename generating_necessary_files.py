import os
import shutil
import sys
import pandas as pd
import pickle
from scripts.topology import (
    ensure_dir,
    get_transcripts_and_sequences,
    align_protein_sequences,
    create_membrane_topology_objects
)
from scripts.expression import getting_expression_data, create_expression_figure_objects
import shelve
import requests
import subprocess
import anndata as ad
from DeepLoc2.deeploc2 import predict
import time 

deeptmhmm_path = os.path.abspath("./DeepTMHMM")
sys.path.append(deeptmhmm_path)

# pyrefly: ignore [missing-import]
from DeepTMHMM.predict_api import load_models, predict_from_fasta

import torch
import DeepLoc2.model
import DeepLoc2.deeploc2
import DeepTMHMM.predict_api
import pytorch_lightning as pl
import torch.nn as nn
import pkg_resources
from esm import pretrained
from DeepTMHMM.utils import hash_aa_string

# Global variables for monkeypatching / caching
_cached_esm_model = None
_cached_deeploc_model = None
_cached_embeddings = {}
_current_sequence = None
model_loading_lock = None

# 1. Patch DeepTMHMM's generate_esm_embeddings to store computed embeddings in memory
original_generate_esm_embeddings = DeepTMHMM.predict_api.generate_esm_embeddings

def patched_generate_esm_embeddings(sequences, esm_embeddings_dir, esm_model, esm_alphabet, repr_layers=33, chunk_size=3):
    original_generate_esm_embeddings(sequences, esm_embeddings_dir, esm_model, esm_alphabet, repr_layers, chunk_size)
    global _cached_embeddings
    for seq in sequences:
        h = hash_aa_string(seq)
        path = os.path.join(esm_embeddings_dir, h)
        if os.path.exists(path):
            _cached_embeddings[seq] = torch.load(path)

DeepTMHMM.predict_api.generate_esm_embeddings = patched_generate_esm_embeddings

# 2. Patch DeepLoc2's run_model_esm1b to pass current sequence and reuse instantiated model
def patched_run_model_esm1b(embed_dataloader, args, test_df):
    global _current_sequence, _cached_deeploc_model
    multilabel_dict = {}
    signaltype_dict = {}
    attn_dict = {}
    memtype_dict = {}
    with torch.no_grad():
        if _cached_deeploc_model is None:
            _cached_deeploc_model = DeepLoc2.model.ESM1bE2E().to(args.device)
        model = _cached_deeploc_model
        for i, (toks, lengths, np_mask, labels) in enumerate(embed_dataloader):
              seq_id = labels[0]
              _current_sequence = test_df.loc[test_df['ACC'] == seq_id, 'Sequence'].values[0]
              ml_out, attn_out, st_out, mt_out = model(toks, lengths, np_mask)
              multilabel_dict[labels[0]] = ml_out
              signaltype_dict[labels[0]] = st_out
              attn_dict[labels[0]] = attn_out
              memtype_dict[labels[0]] = mt_out

    multilabel_df = pd.DataFrame(multilabel_dict.items(), columns=['ACC', 'multilabel'])
    signaltype_df = pd.DataFrame(signaltype_dict.items(), columns=['ACC', 'signaltype'])
    attn_df = pd.DataFrame(attn_dict.items(), columns=['ACC', 'Attention'])
    memtype_df = pd.DataFrame(memtype_dict.items(), columns=['ACC', 'memtype'])

    pred_df = test_df.merge(multilabel_df).merge(signaltype_df).merge(attn_df).merge(memtype_df)
    return pred_df

DeepLoc2.deeploc2.run_model_esm1b = patched_run_model_esm1b

# 3. Patch ESM1bE2E init to reuse pre-loaded ESM model
def patched_ESM1bE2E_init(self):
    pl.LightningModule.__init__(self)
    global _cached_esm_model
    if _cached_esm_model is None:
        
        _cached_esm_model, _ = pretrained.load_model_and_alphabet("esm1b_t33_650M_UR50S")
    self.embedding_func = _cached_esm_model.eval()
    self.subcel_clfs = nn.ModuleList([
        DeepLoc2.model.ESM1bFrozen.load_from_checkpoint(
            pkg_resources.resource_filename("DeepLoc2", f"models/models_esm1b/{i}_1Layer.ckpt"),
            map_location="cpu"
        ).eval() for i in range(5)
    ])
    self.signaltype_clfs = nn.ModuleList([
        DeepLoc2.model.SignalTypeMLP.load_from_checkpoint(
            pkg_resources.resource_filename("DeepLoc2", f"models/models_esm1b/signaltype/{i}.ckpt"),
            map_location="cpu"
        ).eval() for i in range(5)
    ])
    self.memtype_clfs = nn.ModuleList([
        DeepLoc2.model.ESM1B_memtype.load_from_checkpoint(
            pkg_resources.resource_filename("DeepLoc2", f"models/models_esm1b/{i}_1Layer_memtype.ckpt"),
            map_location="cpu"
        ).eval() for i in range(5)
    ])

DeepLoc2.model.ESM1bE2E.__init__ = patched_ESM1bE2E_init
if hasattr(DeepLoc2.deeploc2, 'ESM1bE2E'):
    DeepLoc2.deeploc2.ESM1bE2E.__init__ = patched_ESM1bE2E_init

# 4. Patch ESM1bE2E forward to use cached embedding if available
def patched_ESM1bE2E_forward(self, toks, lens, non_mask):
    global _current_sequence, _cached_embeddings
    device = self.device
    if _current_sequence in _cached_embeddings:
        x = _cached_embeddings[_current_sequence].to(device).unsqueeze(0).float()
    else:
        x = self.embedding_func(toks.to(device), repr_layers=[33])["representations"][33][:, 1:-1].float()
    x_loc_preds, x_signal_preds, x_attnss, x_memtype_preds = [], [], [], []
    for i in range(5):
      x_pred, x_pool, x_attns = self.subcel_clfs[i].predict(x, lens.to(device), non_mask[:, 1:-1].to(device))
      x_loc_preds.append(torch.sigmoid(x_pred))
      x_attnss.append(x_attns)
      x_signal_pred = torch.sigmoid(self.signaltype_clfs[i](torch.cat((x_pool, torch.sigmoid(x_pred)), dim=1)))
      x_signal_preds.append(x_signal_pred)
      x_memtype_pred, _ = self.memtype_clfs[i].forward(x, lens.to(device), non_mask[:, 1:-1].to(device))
      x_memtype_preds.append(torch.sigmoid(x_memtype_pred))
    return torch.stack(x_loc_preds).mean(0).cpu().numpy(), torch.stack(x_attnss).mean(0).cpu().numpy(), torch.stack(x_signal_preds).mean(0).cpu().numpy(), torch.stack(x_memtype_preds).mean(0).cpu().numpy()

DeepLoc2.model.ESM1bE2E.forward = patched_ESM1bE2E_forward
if hasattr(DeepLoc2.deeploc2, 'ESM1bE2E'):
    DeepLoc2.deeploc2.ESM1bE2E.forward = patched_ESM1bE2E_forward

_worker_stats = None

def process_gene_worker(args):
    global _worker_stats, _cached_esm_model, model_loading_lock
    gene_id, gene_names, out_dir_base = args
    worker_out_dir = os.path.join(out_dir_base, f"worker_{gene_id}")
    ensure_dir(worker_out_dir)
    
    # Restrict PyTorch to single-threaded execution per worker to avoid core over-subscription and hybrid CPU bottlenecks
    import torch
    torch.set_num_threads(1)
    
    # Safely load models sequentially using a multiprocessing Lock to prevent system lockup
    if _worker_stats is None:
        if model_loading_lock is not None:
            model_loading_lock.acquire()
        try:
            if _worker_stats is None:
                print(f"[{gene_id}] Worker process initializing and loading models...")
                _worker_stats = load_models()
                _cached_esm_model = _worker_stats[0]
        finally:
            if model_loading_lock is not None:
                model_loading_lock.release()
    
    t_start = time.time()
    
    try:
        print(f"Worker started for {gene_id}")

        # Getting transcripts and sequences
        transcripts_id, mapping, one_isoform = get_transcripts_and_sequences(gene_id, worker_out_dir, session)
        
        if transcripts_id is None:
            return gene_id, gene_names, None, None, None, None

        # Aligning protein sequences
        if not one_isoform:
            base_dir = os.path.dirname(os.path.abspath(__file__))
            align_protein_sequences(mapping, worker_out_dir, base_dir)
        else:
            sequences = open(f"{worker_out_dir}/isoforms.fasta", "r").read()
            with open(f"{worker_out_dir}/aligned_sequences.fasta", "w") as file:
                for i, line in enumerate(sequences.split("\n")):
                    if i == 0: file.write(">" + mapping[line[1:]] + "\n")
                    elif line == "": continue
                    elif line[0] == ">": file.write("\n>" + mapping[line[1:]] + "\n")
                    else: file.write(line.strip())

        # Clear deeptmhmm output
        if os.path.isdir(f"{worker_out_dir}/DeepTMHMM_results"):   
            shutil.rmtree(f"{worker_out_dir}/DeepTMHMM_results")

        # Running DeepTMHMM
        t = time.time()
        predict_from_fasta(f"{worker_out_dir}/isoforms.fasta", f"{worker_out_dir}/DeepTMHMM_results", _worker_stats)
        print(f"DeepTMHMM: {time.time() - t}")

        # Creating membrane topology objects
        membrane_topology_object = create_membrane_topology_objects(transcripts_id, mapping, worker_out_dir)

        # Getting expression data
        expression_normalized_df = getting_expression_data(tpm_expression_df, metadata_df, transcripts_id, mapping)
        if expression_normalized_df is None:
            TCGA_GTEx_plotting_data = None
        else:
            TCGA_GTEx_plotting_data = create_expression_figure_objects(expression_normalized_df)

        # Running deeploc2
        import sys
        old_argv = sys.argv
        sys.argv = [
            "deeploc2", 
            "-f", f"{worker_out_dir}/isoforms.fasta", 
            "-o", f"{worker_out_dir}/deeploc2_output"
        ]
        predict()
        sys.argv = old_argv

        # Reading and cleaning deeploc2 output
        deeploc2_output_name = os.listdir(f"{worker_out_dir}/deeploc2_output")[0]
        deeploc2_output = pd.read_csv(f"{worker_out_dir}/deeploc2_output/{deeploc2_output_name}", index_col=0)
        
        columns_to_remove = ["Localizations", "Signals", "Membrane types"]
        deeploc2_output = deeploc2_output.drop(columns=columns_to_remove)
        deeploc2_output.index = deeploc2_output.index.map(lambda x: mapping[x])

        return gene_id, gene_names, mapping, membrane_topology_object, TCGA_GTEx_plotting_data, deeploc2_output

    except Exception as e:
        print(f"[{gene_id}] Worker error: {e}")
        return gene_id, gene_names, "ERROR", str(e), None, None

    finally:
        # Clean up temporary folder
        if os.path.exists(worker_out_dir):
            shutil.rmtree(worker_out_dir)
        # Clear in-memory cache to prevent memory leaks across runs
        _cached_embeddings.clear()

## Variables

# Output directory
out_dir = "./files"
out_dir_for_plots = "./files_for_plots"

ensure_dir(out_dir_for_plots)

# Proteins we want to show
genes = pd.read_csv("./scripts/files_genes/genes.txt", index_col=0)
genes_id = genes.index.unique()
with open("./scripts/files_genes/already_explored_genes.pkl", "rb") as f:
    already_explored_genes = pickle.load(f)

# Token for biolib
os.environ["BIOLIB_TOKEN"] = "nqCGENf1OAEX405liWxueG9qV9cv6mSETAU3OO8DkjY"

# Create a persistent session outside the function to reuse connections
email = "s242830@dtu.dk"
session = requests.Session()
session.headers.update({"User-Agent": f"PythonBioScript/1.0 ({email})"})

# Batch pre-fetch unexplored sequences from Ensembl
print(f"Getting all unexplored genes' names...")
already_explored_set = set(already_explored_genes)
gene_info = {}
for row in genes.itertuples():
    gid = row[0]
    name = row[1]
    synonym = row[2]
    if gid not in gene_info:
        gene_info[gid] = ([], [])
    gene_info[gid][0].append(name)
    gene_info[gid][1].append(synonym)

gid_to_names = {}
for gid, (names, synonyms) in gene_info.items():
    gene_names = gid
    if len(names) > 1:
        first_name = names[0]
        if pd.notna(first_name):
            gene_names += "|" + str(first_name)
        for syn in synonyms:
            if pd.notna(syn):
                gene_names += "|" + str(syn)
    else:
        name = names[0]
        syn = synonyms[0]
        if pd.notna(name):
            gene_names += "|" + str(name)
            if pd.notna(syn):
                gene_names += "|" + str(syn)
    gid_to_names[gid] = gene_names

print(f"Getting unexplored genes...")
unexplored_genes = []
for gid in genes_id:
    gene_names = gid_to_names.get(gid, gid)
    if gene_names not in already_explored_set:
        unexplored_genes.append((gid, gene_names))

# Loading expression datasets
print("Loading expression datasets...")
tpm_expression_df = ad.read_h5ad(
    "./scripts/files_genes/transcript_counts_sparse.h5ad", backed="r"
)
metadata_df = pd.read_csv(
    "./scripts/files_genes/TcgaTargetGTEX_phenotype.txt",
    index_col=0,
    sep="\t",
    usecols=["sample", "_study", "detailed_category"],
    encoding="latin1",
)

# Loading DeepTMHMM models (will be loaded on-demand in workers)
print("Loading DeepTMHMM models (on-demand in workers)...")

# Running the analysis in parallel using multiprocessing
if __name__ == "__main__":
    import multiprocessing
    
    # unexplored_genes contains tuples of (gene_id, gene_names)
    tasks = [(gene_id, gene_names, out_dir) for gene_id, gene_names in unexplored_genes]
        
    # Explicitly use 'fork' start method to inherit preloaded variables (stats, expression data, etc.)
    mp_ctx = multiprocessing.get_context("fork")
    
    model_loading_lock = mp_ctx.Lock()
    
    with mp_ctx.Pool(processes=50) as pool:
        # Use imap_unordered to process tasks as they finish, writing databases and pickle file immediately
        for result in pool.imap_unordered(process_gene_worker, tasks):
            gene_id, gene_names, mapping, topology, expression, deeploc = result
            
            if mapping == "ERROR":
                print(f"[{gene_id}] Worker failed with error: {topology}. Skipping and writing None to databases.")
                
                # Update the explored genes list and save it
                already_explored_genes.append(gene_names)
                
                print(f"[{gene_id}] {len(already_explored_genes)} / {len(genes_id)} | {round(len(already_explored_genes) / len(genes_id) * 100,3)}%")
                
                # Write None to all 4 databases to support graceful front-end displaying
                with shelve.open(out_dir_for_plots + "/transcripts_to_isoforms_mapping") as db:
                    db[gene_names] = None
                with shelve.open(out_dir_for_plots + "/membrane_topology_objects") as db:
                    db[gene_names] = None
                with shelve.open(out_dir_for_plots + "/TCGA_GTEx_plotting_data") as db:
                    db[gene_names] = None
                with shelve.open(out_dir_for_plots + "/deeploc2_output") as db:
                    db[gene_names] = None
                    
                # Save the already explored genes pickle file immediately
                with open("./scripts/files_genes/already_explored_genes.pkl", "wb") as f:
                    pickle.dump(already_explored_genes, f)
                continue
                
            # Update the explored genes list and save it
            already_explored_genes.append(gene_names)
            
            print(f"[{gene_id}] {len(already_explored_genes)} / {len(genes_id)} | {round(len(already_explored_genes) / len(genes_id) * 100,3)}%")
            
            if mapping is None:
                # Gene encodes no valid protein
                with shelve.open(out_dir_for_plots + "/transcripts_to_isoforms_mapping") as db:
                    db[gene_names] = None
                with shelve.open(out_dir_for_plots + "/membrane_topology_objects") as db:
                    db[gene_names] = None
                with shelve.open(out_dir_for_plots + "/TCGA_GTEx_plotting_data") as db:
                    db[gene_names] = None
                with shelve.open(out_dir_for_plots + "/deeploc2_output") as db:
                    db[gene_names] = None
            else:
                with shelve.open(out_dir_for_plots + "/transcripts_to_isoforms_mapping") as db:
                    db[gene_names] = mapping
                with shelve.open(out_dir_for_plots + "/membrane_topology_objects") as db:
                    db[gene_names] = topology
                with shelve.open(out_dir_for_plots + "/TCGA_GTEx_plotting_data") as db:
                    db[gene_names] = expression
                with shelve.open(out_dir_for_plots + "/deeploc2_output") as db:
                    db[gene_names] = deeploc
                    
            # Save the already explored genes pickle file immediately
            with open("./scripts/files_genes/already_explored_genes.pkl", "wb") as f:
                pickle.dump(already_explored_genes, f)