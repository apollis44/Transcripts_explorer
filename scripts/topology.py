import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import seaborn as sns
import numpy as np
import os
import json
import biolib
from scripts.extract_sequences import (
    fetch_transcripts, 
    fetch_protein_sequence, 
    align_sequences
)

def ensure_dir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def get_transcripts_and_sequences(ensembl_id, output_dir):
    ensure_dir(output_dir)
    
    # 1. Fetch Transcripts
    transcripts = fetch_transcripts(ensembl_id)
    transcripts_ids = []
    for transcript in transcripts:
        transcripts_ids.append(transcript['id'])

    transcripts_ids.sort()

    # 2. Fetch Protein Sequences
    protein_sequences = fetch_protein_sequence(transcripts_ids)

    unique_sequences = list(dict.fromkeys(protein_sequences.values()))
    transcripts_mapping = {}
    with open(f"{output_dir}/isoforms.fasta", "w") as fasta_file:
        for i, seq in enumerate(unique_sequences):
            ids = [k for k, v in protein_sequences.items() if v == seq]
            header = ">"+ "Isoform_" + str(i + 1)
            fasta_file.write(f"{header}\n{seq}\n")
            for id in ids:
                transcripts_mapping[id] = header.split(">")[-1]
            
    return transcripts_ids, transcripts_mapping

def align_protein_sequences(email, output_dir):
    ensure_dir(output_dir)
    
    # Align isoforms
    sequences = open(f"{output_dir}/isoforms.fasta", "r").read()
    alignment = align_sequences(sequences, email)
    with open(f"{output_dir}/aligned_sequences.fasta", "w") as file:
        for i, line in enumerate(alignment.split("\n")):
            if i == 0: file.write(line + "\n")
            elif line == "": continue
            elif line[0] == ">": file.write("\n" + line + "\n")
            else: file.write(line.strip())
            
    return

def run_deeptmhmm(output_dir):
    ensure_dir(f"{output_dir}/DeepTMHMM_results/")

    deeptmhmm = biolib.load('DTU/DeepTMHMM')

    print("Running DeepTMHMM...")
    
    # Check if already run and if it has, empty the directory
    if os.path.exists(f"{output_dir}/DeepTMHMM_results/"):
        files_to_delete = os.listdir(f"{output_dir}/DeepTMHMM_results")
        for file in files_to_delete:
            os.remove(os.path.join(f"{output_dir}/DeepTMHMM_results", file))    
        
    try:                 
        job = deeptmhmm.cli(
            args=f"--fasta {output_dir}/isoforms.fasta", 
        )
        job.save_files(f"{output_dir}/DeepTMHMM_results/")
    except Exception as e:
        print(f"Error running DeepTMHMM: {e}")

    print("DeepTMHMM run completed.")
    return

def create_membrane_topology_objects(mapping, output_dir):
            
    membrane_topology_file = open(f"{output_dir}/DeepTMHMM_results/predicted_topologies.3line").readlines()
    membrane_topology = pd.DataFrame(index=list(set(mapping.values())), columns=["sequence", "topology"])

    # Extracting the alignment and adding it to the dataframe
    aligned_sequences = open(f"{output_dir}/aligned_sequences.fasta").readlines()
    for i, line in enumerate(aligned_sequences):
        if i % 2 == 1:
            isoform_id = aligned_sequences[i - 1].replace(">", "").strip()
            membrane_topology.at[isoform_id, "sequence"] = line

    # Extracting the topology and adding it to the dataframe
    for i, line in enumerate(membrane_topology_file):
        full_topology = ""
        if i % 3 == 0:
            isoform_id = line.split(" ")[0].replace(">", "").strip()
            
        if i % 3 == 2: # The 3-line file is formatted as: [sequence name] [sequence] [topology]
            # Add the topology matching to the alignment (- will be matched with -)
            topology = line.strip()
            aligned_seq = membrane_topology.at[isoform_id, "sequence"].strip()
            j = 0
            for char in aligned_seq:
                if char == "-":
                    full_topology += "-"
                else:
                    full_topology += topology[j]
                    j += 1
            membrane_topology.at[isoform_id, "topology"] = full_topology

    # Create the data for each sequence
    # The data is a list of (start, width) tuples for each feature (fx. [{'-': [(0, 95), (194, 694)],'E': [(95, 99)])
    sequences_data = []
    for i in range(len(membrane_topology)):
        isoform_id = "Isoform_" + str(i + 1)
        topology = membrane_topology.at[isoform_id, "topology"]
        seq_data = {}
        # Identify features in the topology
        current_feature = None
        start = None
        for pos, char in enumerate(topology):
            if char != current_feature:
                if current_feature is not None:
                    width = pos - start
                    if current_feature not in seq_data:
                        seq_data[current_feature] = []
                    seq_data[current_feature].append((start, width))
                current_feature = char
                start = pos
        if current_feature is not None:
            width = len(topology) - start
            if current_feature not in seq_data:
                seq_data[current_feature] = []
            seq_data[current_feature].append((start, width))
        sequences_data.append(seq_data)
    
    # Save membrane_topology and sequences_data
    membrane_topology.to_csv(output_dir + "/membrane_topology.csv", index=True)
    with open(f"{output_dir}/sequences_data.json", "w") as f:
        json.dump(sequences_data, f, indent=4)
    print(f"Saved sequences_data.json to {output_dir}")

    return sequences_data