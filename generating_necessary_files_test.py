import os
import shutil
import sys
import matplotlib.pyplot as plt
import pandas as pd
import pickle
from scripts.topology import (
    ensure_dir,
    get_transcripts_and_sequences,
    align_protein_sequences,
    run_deeptmhmm,
    create_membrane_topology_objects
)
from scripts.expression import (
    getting_expression_data,
    create_expression_figure_objects
)
import shelve
import requests
import subprocess
import time
import xenaPython as xena

deeptmhmm_path = os.path.abspath("./DeepTMHMM")
sys.path.append(deeptmhmm_path)

from DeepTMHMM.predict_api import (
    load_models,
    predict_from_fasta
)

## Variables

# Output directory
out_dir = "./files"
out_dir_for_plots = "./app/files_for_plots"

# Proteins we want to show
genes = pd.read_csv("./scripts/files_genes/genes.txt", index_col=0)
genes_id = genes.index.unique()
with open("./scripts/files_genes/already_explored_genes.pkl", "rb") as f:
    already_explored_genes = pickle.load(f)

# Token for biolib
os.environ['BIOLIB_TOKEN'] = "nqCGENf1OAEX405liWxueG9qV9cv6mSETAU3OO8DkjY"

# Create a persistent session outside the function to reuse connections
email = "s242830@dtu.dk"
session = requests.Session()
session.headers.update({
    "User-Agent": f"PythonBioScript/1.0 ({email})"
})

# Precomputed heavy xenaPython steps

# Getting the right dataset from Xena
print("Getting the right dataset from Xena...")
host = "https://toil.xenahubs.net" # Public hub with both TCGA and GTEx data
cohort = "TCGA TARGET GTEx" # Cohort name
samples = xena.cohort_samples(host, cohort, None)
dataset = "TcgaTargetGtex_rsem_isoform_tpm" # Dataset name
metadata_dataset = "TcgaTargetGTEX_phenotype.txt"
fields = ["_study", "detailed_category"] # Fields to extract

# We extract the list of all transcripts ids available in the dataset
available_transcripts = xena.dataset_field(host, dataset)

# The dataset contains codes for categorical fields, we convert them to their actual values
# To do so, we fetch the codes from Xena and create a mapping dictionary
codes = xena.field_codes(host, metadata_dataset, fields)

# Loading DeepTMHMM models
print("Loading DeepTMHMM models...")
stats = load_models()

genes_id = ["ENSG00000156738"] # Testing on one gene only (CD20)
# Running the analysis for each protein
for gene_id in genes_id:
    print(gene_id)
    t = time.time()
    gene_names = gene_id

    # Finding gene name and synonyms
    if isinstance(genes.loc[gene_id,"Gene name"], pd.Series):
        gene_names += "|" + genes.loc[gene_id,"Gene name"].iloc[0]

        for synonym in genes.loc[gene_id,"Gene Synonym"].dropna(axis=0).to_list():
            gene_names += "|" + synonym

    elif not pd.isna(genes.loc[gene_id,"Gene name"]):
        gene_names += "|" + genes.loc[gene_id,"Gene name"]
        
        if not pd.isna(genes.loc[gene_id,"Gene Synonym"]):
            gene_names += "|" + genes.loc[gene_id,"Gene Synonym"]

    already_explored_genes.append(gene_names)

    print(f"Processing {gene_id}...")
    ensure_dir(out_dir)

    print("Time to get the gene names: " + str(time.time() - t))
    t = time.time()

    # Getting transcripts and sequences
    transcripts_id, mapping, transcripts_length, one_isoform = get_transcripts_and_sequences(gene_id, out_dir, session)
    
    if transcripts_id == None: # The gene encodes no valid protein
        with open(out_dir + "/already_explored_genes.pkl", "wb") as f:
            pickle.dump(already_explored_genes, f)
        continue

    print("Time to get the transcripts and sequences: " + str(time.time() - t))
    t = time.time()

    # Aligning protein sequences
    if not one_isoform:
        base_dir = os.path.dirname(os.path.abspath(__file__))
        align_protein_sequences(mapping, out_dir, base_dir)
    else:
        sequences = open(f"{out_dir}/isoforms.fasta", "r").read()
        with open(f"{out_dir}/aligned_sequences.fasta", "w") as file:
            for i, line in enumerate(sequences.split("\n")):
                if i == 0: file.write(">" + mapping[line[1:]] + "\n")
                elif line == "": continue
                elif line[0] == ">": file.write("\n>" + mapping[line[1:]] + "\n")
                else: file.write(line.strip())

    print("Time to align the protein sequences: " + str(time.time() - t))
    t = time.time()

    # Clear deeptmhmm output
    if os.path.isdir(f"{out_dir}/DeepTMHMM_results"):   
        shutil.rmtree(f"{out_dir}/DeepTMHMM_results")

    # Running DeepTMHMM
    predict_from_fasta(f"{out_dir}/isoforms.fasta", f"{out_dir}/DeepTMHMM_results", stats)

    print("Time to run DeepTMHMM: " + str(time.time() - t))
    t = time.time()

    # Storing transcripts to isoforms mapping
    with shelve.open(out_dir_for_plots + "/transcripts_to_isoforms_mapping") as db:
        db[gene_names] = mapping

    # Creating membrane topology objects
    membrane_topology_object = create_membrane_topology_objects(transcripts_id, mapping, out_dir)

    print("Time to create the membrane topology objects: " + str(time.time() - t))
    t = time.time()

    # Storing membrane topology objects
    with shelve.open(out_dir_for_plots + "/membrane_topology_objects") as db:
        db[gene_names] = membrane_topology_object

    # Getting expression data
    expression_normalized_df = getting_expression_data(transcripts_id, transcripts_length, mapping, samples, available_transcripts, codes)
    if expression_normalized_df is None:
        TCGA_GTEx_plotting_data = None
    else:
        TCGA_GTEx_plotting_data = create_expression_figure_objects(expression_normalized_df)

    print("Time to get the expression data: " + str(time.time() - t))
    t = time.time()

    # Running deeploc2 in command lines
    command = [
        "deeploc2", 
        "-f", f"{out_dir}/isoforms.fasta", 
        "-o", f"{out_dir}/deeploc2_output"
    ]
    subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # Reading deeploc2 output
    deeploc2_output_name = os.listdir(f"{out_dir}/deeploc2_output")[0]
    deeploc2_output = pd.read_csv(f"{out_dir}/deeploc2_output/{deeploc2_output_name}", index_col=0)
    subprocess.run(f"rm -r {out_dir}/deeploc2_output")

    # Cleaning deeploc2 output
    columns_to_remove = ["Localizations", "Signals", "Membrane types"]
    deeploc2_output = deeploc2_output.drop(columns=columns_to_remove)

    print("Time to run deeploc2: " + str(time.time() - t))
    t = time.time()