import os
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

## Variables

# Output directory
out_dir = "./files"
out_dir_for_plots = "./app/files_for_plots"

# Proteins we want to show
genes = pd.read_csv(out_dir + "/genes.txt", index_col=0)
genes_id = genes.index.unique()
with open(out_dir + "/already_explored_genes.pkl", "rb") as f:
    already_explored_genes = pickle.load(f)

# Token for biolib
os.environ['BIOLIB_TOKEN'] = "nqCGENf1OAEX405liWxueG9qV9cv6mSETAU3OO8DkjY"

# Create a persistent session outside the function to reuse connections
email = "s242830@dtu.dk"
session = requests.Session()
session.headers.update({
    "User-Agent": f"PythonBioScript/1.0 ({email})"
})

# Running the analysis for each protein
for gene_id in genes_id:
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

    # We skip genes that were already explored
    if gene_names in already_explored_genes:
        continue

    already_explored_genes.append(gene_names)

    print(f"Processing {gene_id}... ({len(already_explored_genes)} / {len(genes_id)})")
    ensure_dir(out_dir)

    # Getting transcripts and sequences
    transcripts_id, mapping, transcripts_length, one_isoform = get_transcripts_and_sequences(gene_id, out_dir, session)
    
    if transcripts_id == None: # The gene encodes no valid protein
        with open(out_dir + "/already_explored_genes.pkl", "wb") as f:
            pickle.dump(already_explored_genes, f)
        continue

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

    # Running DeepTMHMM
    run_deeptmhmm(out_dir)

    # Storing transcripts to isoforms mapping
    with shelve.open(out_dir_for_plots + "/transcripts_to_isoforms_mapping") as db:
        db[gene_names] = mapping

    # Creating membrane topology objects
    membrane_topology_object = create_membrane_topology_objects(transcripts_id, mapping, out_dir)

    # Storing membrane topology objects
    with shelve.open(out_dir_for_plots + "/membrane_topology_objects") as db:
        db[gene_names] = membrane_topology_object

    # Getting expression data
    expression_normalized_df = getting_expression_data(transcripts_id, transcripts_length, mapping)
    if expression_normalized_df is None:
        TCGA_GTEx_plotting_data = None
    else:
        TCGA_GTEx_plotting_data = create_expression_figure_objects(expression_normalized_df)

    # Storing expression data
    with shelve.open(out_dir_for_plots + "/TCGA_GTEx_plotting_data") as db:
        db[gene_names] = TCGA_GTEx_plotting_data

    # Running deeploc2 in command lines
    subprocess.run(f"deeploc2 -f {out_dir}/isoforms.fasta -o {out_dir}/deeploc2_output", stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # Reading deeploc2 output
    deeploc2_output_name = os.listdir(f"{out_dir}/deeploc2_output")[0]
    deeploc2_output = pd.read_csv(f"{out_dir}/deeploc2_output/{deeploc2_output_name}", index_col=0)
    subprocess.run(f"rm -r {out_dir}/deeploc2_output")

    # Cleaning deeploc2 output
    columns_to_remove = ["Localizations", "Signals", "Membrane types"]
    deeploc2_output = deeploc2_output.drop(columns=columns_to_remove)

    # Storing deeploc2 output
    with shelve.open(out_dir_for_plots + "/deeploc2_output") as db:
        db[gene_names] = deeploc2_output

    # Storing already explored genes
    with open(out_dir + "/already_explored_genes.pkl", "wb") as f:
        pickle.dump(already_explored_genes, f)