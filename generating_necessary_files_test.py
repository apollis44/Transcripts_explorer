import os
import shutil
import sys
import anndata as ad
import pandas as pd
import shelve
from scripts.topology import (
    ensure_dir,
    get_transcripts_and_sequences,
    align_protein_sequences,
    create_membrane_topology_objects,
)
from scripts.expression import (
    getting_expression_data,
    create_expression_figure_objects
)
import requests
import subprocess
import time
from DeepLoc2.deeploc2 import predict

deeptmhmm_path = os.path.abspath("./DeepTMHMM")
sys.path.append(deeptmhmm_path)

from DeepTMHMM.predict_api import (
    load_models,
    predict_from_fasta
)

import torch
import DeepLoc2.model
import DeepLoc2.deeploc2
import DeepTMHMM.predict_api
import pytorch_lightning as pl
import torch.nn as nn
import pkg_resources

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
    from DeepTMHMM.utils import hash_aa_string
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
        from esm import pretrained
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
    print(f"[{gene_id}] Worker started.")
    
    try:
        # Getting transcripts and sequences
        t = time.time()
        transcripts_id, mapping, one_isoform = get_transcripts_and_sequences(gene_id, worker_out_dir, session)
        print(f"[{gene_id}] Time to get the transcripts and sequences: {time.time() - t:.4f}")
        
        if transcripts_id is None:
            print(f"[{gene_id}] No valid protein encoded.")
            return gene_id, gene_names, None, None, None, None

        # Aligning protein sequences
        t = time.time()
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
        print(f"[{gene_id}] Time to align the protein sequences: {time.time() - t:.4f}")

        # Clear deeptmhmm output
        if os.path.isdir(f"{worker_out_dir}/DeepTMHMM_results"):   
            shutil.rmtree(f"{worker_out_dir}/DeepTMHMM_results")

        # Running DeepTMHMM
        t = time.time()
        predict_from_fasta(f"{worker_out_dir}/isoforms.fasta", f"{worker_out_dir}/DeepTMHMM_results", _worker_stats)
        print(f"[{gene_id}] Time to run DeepTMHMM: {time.time() - t:.4f}")

        # Creating membrane topology objects
        t = time.time()
        membrane_topology_object = create_membrane_topology_objects(transcripts_id, mapping, worker_out_dir)
        print(f"[{gene_id}] Time to create the membrane topology objects: {time.time() - t:.4f}")

        # Getting expression data
        t = time.time()
        expression_normalized_df = getting_expression_data(tpm_expression_df, metadata_df, transcripts_id, mapping)
        if expression_normalized_df is None:
            TCGA_GTEx_plotting_data = None
        else:
            TCGA_GTEx_plotting_data = create_expression_figure_objects(expression_normalized_df)
        print(f"[{gene_id}] Time to get the expression data: {time.time() - t:.4f}")

        # Running deeploc2
        t = time.time()
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
        print(f"[{gene_id}] Time to run deeploc2: {time.time() - t:.4f}")

        print(f"[{gene_id}] Worker finished in {time.time() - t_start:.2f} seconds.")
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

# Proteins we want to show
genes = pd.read_csv("./scripts/files_genes/genes.txt", index_col=0)
genes_id = genes.index.unique()

# Token for biolib
os.environ['BIOLIB_TOKEN'] = "nqCGENf1OAEX405liWxueG9qV9cv6mSETAU3OO8DkjY"

# Create a persistent session outside the function to reuse connections
email = "s242830@dtu.dk"
session = requests.Session()
session.headers.update({
    "User-Agent": f"PythonBioScript/1.0 ({email})"
})

# Loading expression datasets
print("Loading expression datasets...")
tpm_expression_df = ad.read_h5ad("./scripts/files_genes/transcript_counts_sparse.h5ad", backed="r")
metadata_df = pd.read_csv(f"./scripts/files_genes/TcgaTargetGTEX_phenotype.txt", 
                          index_col=0, 
                          sep="\t", 
                          usecols=["sample", "_study", "detailed_category"], 
                          encoding="latin1" 
)

# Loading DeepTMHMM models (will be loaded on-demand in workers)
print("Loading DeepTMHMM models (on-demand in workers)...")

genes_id = ["ENSG00000177455", "ENSG00000156738", "ENSG00000003056", "ENSG00000261857"]
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

# Running the analysis in parallel using multiprocessing
if __name__ == "__main__":
    import multiprocessing
    
    # We dynamically select up to 4 genes that are in the cache and not None,
    # fallback to ["ENSG00000156738"] if none are found.
    test_genes = ["ENSG00000182670","ENSG00000291738","ENSG00000291791", "ENSG00000285479"]
    # with shelve.open("./.cache_ensembl/sequences_cache") as db:
    #     for gid in genes.index.unique():
    #         if gid in db and db[gid] is not None:
    #             test_genes.append(gid)
    #             if len(test_genes) == 4:
    #                 break
                    
    if not test_genes:
        test_genes = ["ENSG00000156738"]
        
    print(f"Testing multiprocessing on genes: {test_genes}")
    
    tasks = [(gid, gid_to_names.get(gid, gid), out_dir) for gid in test_genes]
    
    print(f"Starting multiprocessing pool with {min(len(tasks), 4)} workers...")
    t_multi = time.time()
    
    # Explicitly use 'fork' start method to inherit preloaded variables (stats, expression data, etc.)
    mp_ctx = multiprocessing.get_context("fork")
    
    model_loading_lock = mp_ctx.Lock()
    
    with mp_ctx.Pool(processes=min(len(tasks), 3)) as pool:
        results = pool.map(process_gene_worker, tasks)
        
    print(f"\nAll worker processes completed in {time.time() - t_multi:.2f} seconds.\n")
    
    # Print status summary in the parent process without saving to any databases
    for gene_id, gene_names, mapping, topology, expression, deeploc in results:
        if mapping == "ERROR":
            print(f"[{gene_id}] Worker failed with error: {topology}")
        else:
            print(f"[{gene_id}] Completed processing {gene_names} successfully!")
                
    print("Multiprocessing test pipeline execution completed successfully!")