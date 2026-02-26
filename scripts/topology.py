import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import seaborn as sns
import numpy as np
import os
import json
import biolib
from .extract_sequences import (
    fetch_transcripts, 
    fetch_protein_sequence, 
    fetch_cdna_length
)
import time
import shutil
import subprocess
from io import StringIO
from Bio import SeqIO

def ensure_dir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def get_transcripts_and_sequences(ensembl_id, output_dir, session):
    ensure_dir(output_dir)
    
    # 1. Fetch Transcripts
    transcripts = fetch_transcripts(ensembl_id, session)

    if transcripts == []:
        return None, None, None, None

    transcripts_ids = []
    for transcript in transcripts:
        transcripts_ids.append(transcript['id'])

    transcripts_ids.sort()

    # 2. Fetch Transcripts Lengths
    transcripts_lengths = fetch_cdna_length(transcripts_ids, session)

    # 3. Fetch Protein Sequences
    protein_sequences = fetch_protein_sequence(transcripts_ids, session)

    unique_sequences = list(dict.fromkeys(protein_sequences.values()))
    if len(unique_sequences) == 1:
        one_isoform = True
    else:
        one_isoform = False

    transcripts_mapping = {}
    with open(f"{output_dir}/isoforms.fasta", "w") as fasta_file:
        for i, seq in enumerate(unique_sequences):
            ids = [k for k, v in protein_sequences.items() if v == seq]
            header = ">" + ids[0]
            fasta_file.write(f"{header}\n{seq}\n")
            for id in ids:
                transcripts_mapping[id] = "<br>".join(ids)
            
    return transcripts_ids, transcripts_mapping, transcripts_lengths, one_isoform

def align_protein_sequences(mapping, output_dir, base_dir):
    ensure_dir(output_dir)
    
    # Align isoforms
    fasta_input = f"{output_dir}/isoforms.fasta"
    fasta_output = f"{output_dir}/aligned_sequences.fasta"

    # Try to find 'mafft' in the system PATH (works for Linux/Mac/Conda)
    mafft_executable = shutil.which("mafft")

    # If it's not in the system PATH, use your local Windows version
    if not mafft_executable:
        base_dir = os.path.dirname(os.path.abspath(__file__))
        mafft_executable = os.path.join(base_dir, "..", "mafft-win", "mafft.bat")

    mafft_args = [
        mafft_executable, 
        "--localpair", 
        "--maxiterate", "1000", 
        "--quiet", 
        fasta_input
    ]

    # Run MAFFT and capture stdout in memory (text=True decodes it to a string)
    result = subprocess.run(mafft_args, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, check=True)

    records = list(SeqIO.parse(StringIO(result.stdout), "fasta"))

    # Update the IDs using your mapping dictionary
    for record in records:
        # Check if the ID exists in your mapping to avoid KeyError
        if record.id in mapping:
            record.id = mapping[record.id]
            # Crucial step: clear the description so Biopython doesn't 
            # accidentally append the old header info next to your new ID!
            record.description = "" 

    # Write them to your file without line breaks (fasta-2line)
    output_file = f"{output_dir}/aligned_sequences.fasta"
    with open(output_file, "w") as file:
        SeqIO.write(records, file, "fasta-2line")

    return

def run_deeptmhmm(output_dir):
    ensure_dir(f"{output_dir}/DeepTMHMM_results/")

    deeptmhmm = biolib.load('DTU/DeepTMHMM') 

    local_temp_fasta = "isoforms_temp.fasta"
    shutil.copy(f"{output_dir}/isoforms.fasta", local_temp_fasta)

    # Check if already run and if it has, empty the directory
    if os.path.exists(f"{output_dir}/DeepTMHMM_results/"):
        files_to_delete = os.listdir(f"{output_dir}/DeepTMHMM_results")
        for file in files_to_delete:
            os.remove(os.path.join(f"{output_dir}/DeepTMHMM_results", file))    
        
    try:
        job = deeptmhmm.cli(args=f'--fasta {local_temp_fasta}', blocking=False)

        while job.get_status() == 'in_progress':
            time.sleep(1)

        job.save_files(f"{output_dir}/DeepTMHMM_results/")
    except Exception as e:
        print(f"Error running DeepTMHMM: {e}")

    finally:
        os.remove(local_temp_fasta)

    return

def create_membrane_topology_objects(transcripts_id, mapping, output_dir):
            
    membrane_topology_file = open(f"{output_dir}/DeepTMHMM_results/predicted_topologies.3line").readlines()
    membrane_topology = pd.DataFrame(index=transcripts_id, columns=["sequence", "topology"])

    # Extracting the alignment and adding it to the dataframe
    aligned_sequences = open(f"{output_dir}/aligned_sequences.fasta").readlines()
    for i, line in enumerate(aligned_sequences):
        if i % 2 == 1:
            extracted_transcripts_id = aligned_sequences[i - 1][1:].split("<br>")
            extracted_transcripts_id = [transcript_id.strip() for transcript_id in extracted_transcripts_id]
            for transcript_id in extracted_transcripts_id:
                membrane_topology.at[transcript_id, "sequence"] = line.strip()

    # Extracting the topology and adding it to the dataframe
    for i, line in enumerate(membrane_topology_file):
        if i % 3 == 0:
            extracted_transcripts_id = mapping[line.split(" ")[0][1:]].split("<br>")
            extracted_transcripts_id = [transcript_id.strip() for transcript_id in extracted_transcripts_id]
            
        if i % 3 == 2: # The 3-line file is formatted as: [sequence name] [sequence] [topology]
            # Add the topology matching to the alignment (- will be matched with -)
            for transcript_id in extracted_transcripts_id:
                full_topology = ""
                topology = line.strip()
                aligned_seq = membrane_topology.at[transcript_id, "sequence"]
                j = 0
                for char in aligned_seq:
                    if char == "-":
                        full_topology += "-"
                    else:
                        full_topology += topology[j]
                        j += 1
                membrane_topology.at[transcript_id, "topology"] = full_topology

    # We will now sort the dataframe to have identical isoforms be next to each other
    isoforms = list(dict.fromkeys(mapping.values()))
    order = [transcript_id for isoform in isoforms for transcript_id in isoform.split("<br>")]
    membrane_topology = membrane_topology.loc[order]

    # Create the data for each sequence
    # The data is a list of (start, width) tuples for each feature (fx. [{'-': [(0, 95), (194, 694)],'E': [(95, 99)])
    sequences_data = []
    for i in range(len(membrane_topology)):
        transcript_id = membrane_topology.index[i]
        topology = membrane_topology.at[transcript_id, "topology"]
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

    return sequences_data