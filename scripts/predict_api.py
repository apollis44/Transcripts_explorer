import numpy as np
import io
import time
import argparse
from Bio import SeqIO
import torch
import sys
from collections import OrderedDict
import os
from tqdm import tqdm
from esm import Alphabet, FastaBatchedDataset, ProteinBertModel, pretrained
from hashlib import md5
from glob import glob

from utils import type_id_to_string, gff3, chunk_with_constraints, load_model_from_disk, prot_type_to_display_prot_type, write_probabilities_to_file, hash_aa_string
from experiments.tmhmm3.tm_util import original_labels_to_fasta, is_topologies_equal, crf_states, max_alpha_membrane_length, max_beta_membrane_length, max_signal_length

import warnings
from torch.serialization import SourceChangeWarning
warnings.filterwarnings("ignore", category=SourceChangeWarning)

base_dir = os.path.dirname(os.path.realpath(__file__))

def generate_esm_embeddings(sequences, esm_embeddings_dir, esm_model, esm_alphabet, repr_layers=33, chunk_size=3):

    os.makedirs(esm_embeddings_dir, exist_ok=True)

    with torch.no_grad():
        if torch.cuda.is_available():
            torch.cuda.set_device('cuda:0')
            esm_model = esm_model.cuda()
        esm_model.eval()

        batch_converter = esm_alphabet.get_batch_converter()
        for idx, seq in enumerate(sequences):

            # if os.path.isfile(f'{esm_embeddings_dir}/{hash_aa_string(seq)}'):
            #     print("Already processed sequence")
            #     continue

            #print(f"Average sequence length in chunk: {sum( map(len, sequence_chunk) ) / len(sequence_chunk)} \n")

            embedding_seq_time = time.time()

            seqs = list([("seq", s) for s in [seq]])
            labels, strs, toks = batch_converter(seqs)
            repr_layers_list = [
                (i + esm_model.num_layers + 1) % (esm_model.num_layers + 1) for i in range(repr_layers)
            ]

            if torch.cuda.is_available():
                toks = toks.to(device="cuda", non_blocking=True)

            minibatch_max_length = toks.size(1)

            tokens_list = []
            end = 0
            while end <= minibatch_max_length:
                start = end
                end = start + 1022
                if end <= minibatch_max_length:
                    # we are not on the last one, so make this shorter
                    end = end - 300
                tokens = esm_model(toks[:, start:end], repr_layers=repr_layers_list, return_contacts=False)[
                    "representations"][repr_layers - 1]
                tokens_list.append(tokens)

            out = torch.cat(tokens_list, dim=1)

            # set nan to zeros
            out[out != out] = 0.0

            res = out.transpose(0, 1)[1:-1]
            seq_embedding = res[:, 0]
            output_file = open(f'{esm_embeddings_dir}/{hash_aa_string(seq)}', 'wb')
            torch.save(seq_embedding, output_file)
            output_file.close()
            #
            # for seq_idx, seq in enumerate(sequence_chunk):
            #     output_file = open(f'{esm_embeddings_dir}/{hash_aa_string(seq)}', 'wb')
            #     #print(f"Saving file to: '{esm_embeddings_dir}/{hash_aa_string(seq)}'")
            #     seq_embedding = res[:, seq_idx]
            #     torch.save(seq_embedding, output_file)
            #     output_file.close()

            #print(f'Generated embedding for sequence chunk in {time.time() - embedding_seq_time}s')
            #print(f"Per sequence {(time.time() - embedding_seq_time) / chunk_size}s")

    # Manually delete the ESM model so it does not consume memory when running the DeepTMHMM models
    del esm_model

def load_esm_model():
    model_load_time = time.time()
    #model_data = torch.load('esm1b_model.pt', map_location='cpu')
    #esm_model, esm_alphabet = pretrained.load_model_and_alphabet_core(model_data)
    torch_load_time = time.time()
    esm_args = torch.load(f'{base_dir}/esm_model_args.pt')
    esm_alphabet = torch.load(f'{base_dir}/esm_model_alphabet.pt')
    esm_model_state_dict = torch.load(f'{base_dir}/esm_model_state_dict.pt')
    #print(f'It took {time.time() - torch_load_time} to load torch vars')
    esm_model = ProteinBertModel(
        args=esm_args,
        alphabet=esm_alphabet
    )

    state_dict_time = time.time()
    esm_model.load_state_dict(esm_model_state_dict)
    #print(f'Loaded state dict in {time.time() - state_dict_time}s \n')

    #print(f'Loaded ESM model in {time.time() - model_load_time}s \n')
    return esm_model, esm_alphabet

def load_models(model_paths=None):
    if model_paths is None:
        # REMOVE the leading slashes here
        model_paths = [
            'deeptmhmm_cv_0.model',
            'deeptmhmm_cv_1.model',
            'deeptmhmm_cv_2.model',
            'deeptmhmm_cv_3.model',
            'deeptmhmm_cv_4.model',
        ]
    
    esm_model, esm_alphabet = load_esm_model()
    
    models = []
    for path in model_paths:
        # Now os.path.join will correctly produce: /your/base/dir/deeptmhmm_cv_0.model
        full_path = os.path.join(base_dir, path)
        models.append(load_model_from_disk(full_path))

    for model in models:
        if torch.cuda.is_available():
            model.cuda()
            model.use_gpu = True
        else:
            model.use_gpu = False
        model.eval()
        
    return esm_model, esm_alphabet, models

def predict_from_fasta(fasta_path, output_dir, stats):
    esm_model, esm_alphabet, models = stats

    torch.set_num_threads(os.cpu_count())
    torch.set_num_interop_threads(os.cpu_count() + 2)

    if os.path.exists(output_dir):
        print(f"Error: output directory {output_dir} already exists")
        return
    os.makedirs(output_dir)

    input_str = open(fasta_path, "r").read()

    if not input_str.strip().startswith(">") and input_str.strip() != "":
        input_str = "> Unnamed Protein \n" + input_str

    input_str_buffer = io.StringIO(input_str)

    protein_id_to_sequence = {}
    protein_sequences = []
    for idx, record in enumerate(SeqIO.parse(input_str_buffer, "fasta")):
        if record.seq == '':
            print(f"Error: No sequence for ID {record.id} found.")
            return
        protein_id_to_sequence[str(record.id)] = str(record.seq)
        protein_sequences.append(str(record.seq))

    if len(protein_sequences) == 0:
        print("Error: No sequences found in FASTA input")
        return

    protein_sequence_to_predicted_type = {}
    protein_sequence_to_predicted_labels = {}
    protein_sequence_to_predicted_topology = {}
    protein_sequence_to_original_labels = {}

    chunk_size = 1
    types_count = {"TM": 0, "SP+TM": 0, "SP": 0, "GLOB": 0, "BETA": 0}

    input_sequences = list(reversed(sorted(protein_sequences, key=len)))

    embeddings_dir = f'{output_dir}/embeddings'

    # Generate embeddings
    generate_esm_embeddings(
        sequences=input_sequences,
        esm_embeddings_dir=embeddings_dir,
        esm_model=esm_model,
        esm_alphabet=esm_alphabet,
        chunk_size=1
    )

    for model in models:
        model.esm_embeddings_dir = embeddings_dir

    protein_sequence_chunks = chunk_with_constraints(input_sequences, chunk_size)

    for sequence_chunk in protein_sequence_chunks:
        all_models_predicted_labels = {sequence: {} for sequence in sequence_chunk}
        all_models_predicted_types = {sequence: {} for sequence in sequence_chunk}
        all_models_predicted_topologies = {sequence: {} for sequence in sequence_chunk}
        all_models_crf_loss = {sequence: {} for sequence in sequence_chunk}

        with torch.no_grad():
            model_to_use_idx = None
            for model_idx, model in enumerate(models):
                predicted_labels, predicted_types, predicted_topologies, _, emissions, mask, predicted_crf_labels = model(sequence_chunk)
                padded_predicted_crf_labels = torch.nn.utils.rnn.pad_sequence(predicted_crf_labels)
                crf_loss_batch = model.crf_model(emissions=emissions, tags=padded_predicted_crf_labels, mask=mask, reduction='none')

                for seq_idx, sequence in enumerate(sequence_chunk):
                    type_string = type_id_to_string([predicted_types[seq_idx]])[0]
                    labels = original_labels_to_fasta(predicted_labels[seq_idx])

                    all_models_predicted_types[sequence][model_idx] = type_string
                    all_models_predicted_labels[sequence][model_idx] = predicted_labels[seq_idx]
                    all_models_predicted_topologies[sequence][model_idx] = predicted_topologies[seq_idx]
                    all_models_crf_loss[sequence][model_idx] = crf_loss_batch[seq_idx]

            for sequence in sequence_chunk:
                topology_agreements = {model_idx: 0 for model_idx in range(len(models))}
                for model_idx in range(len(models)):
                    for other_model_idx in range(len(models)):
                        if model_idx == other_model_idx:
                            continue
                        model_topology = all_models_predicted_topologies[sequence][model_idx]
                        other_model_topology = all_models_predicted_topologies[sequence][other_model_idx]
                        if is_topologies_equal(model_topology, other_model_topology):
                            topology_agreements[model_idx] += 1

                max_topology_agreements = max(topology_agreements.values())
                tied_models_indexes = [idx for idx, value in topology_agreements.items() if value == max_topology_agreements]

                if len(tied_models_indexes) > 1:
                    model_to_use_idx = max({idx: loss for idx, loss in all_models_crf_loss[sequence].items() if idx in tied_models_indexes}, key=all_models_crf_loss[sequence].get)
                else:
                    model_to_use_idx = tied_models_indexes[0]

                types_count[all_models_predicted_types[sequence][model_to_use_idx]] += 1
                protein_sequence_to_predicted_type[sequence] = all_models_predicted_types[sequence][model_to_use_idx]
                protein_sequence_to_predicted_labels[sequence] = original_labels_to_fasta(all_models_predicted_labels[sequence][model_to_use_idx])
                protein_sequence_to_original_labels[sequence] = all_models_predicted_labels[sequence][model_to_use_idx]
                protein_sequence_to_predicted_topology[sequence] = all_models_predicted_topologies[sequence][model_to_use_idx]


    if len(protein_sequences) == 1:
        model = models[model_to_use_idx]
        model.use_marg_prob = True
        marg_probs, mask, _ = model.get_emissions_for_decoding(protein_sequences)

    predicted_topologies = ''
    for prot_id, seq in protein_id_to_sequence.items():
        predicted_topologies += f">{prot_id} | {protein_sequence_to_predicted_type[seq]}\n{seq}\n{protein_sequence_to_predicted_labels[seq]}\n"

    with open(f"{output_dir}/predicted_topologies.3line", "w") as outfile:
        outfile.write(predicted_topologies)

    markdown_outfile = open(f'{output_dir}/deeptmhmm_results.md' , 'w')
    gff3_output = "##gff-version 3\n"
    region_count = {"TMhelix": 0, "signal": 0, "inside": 0, "periplasm": 0, "outside": 0, 'Beta sheet': 0}

    print_separator = True
    for seq_idx, prot_seq_id in enumerate(protein_id_to_sequence.keys(), 1):
        sequence = protein_id_to_sequence[prot_seq_id]
        prot_seq_topology = protein_sequence_to_predicted_topology[sequence]
        prot_seq_labels = protein_sequence_to_predicted_labels[sequence]

        tmrs = sum(1 for region in prot_seq_topology if region[1] in (4, 5))
        region_strings = []

        for idx, region in enumerate(prot_seq_topology):
            topology_category = region[1].item()
            end = len(prot_seq_labels) if idx == len(prot_seq_topology) - 1 else int(prot_seq_topology[idx+1][0])
                
            region_map = {0: "inside", 1: "outside", 2: "periplasm", 3: "signal", 4: "TMhelix", 5: "Beta sheet"}
            region_str = region_map.get(topology_category, "unknown region")
            if region_str == "unknown region":
                print("Error: unknown region", file=markdown_outfile)
                
            region_strings.append([region_str, str(int(region[0]) + 1), str(end)])
            if region_str in region_count:
                region_count[region_str] += 1

        if seq_idx == len(protein_sequences):
            print_separator = False

        gff3_output += gff3(prot_seq_id, len(prot_seq_labels), tmrs, region_strings, print_separator)

    with open(f"{output_dir}/TMRs.gff3", "w") as outfile:
        outfile.write(gff3_output)