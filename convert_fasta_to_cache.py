#!/usr/bin/env python3
import os
import shelve
import sys

def parse_fasta(fasta_path):
    """
    Parses the Ensembl peptide FASTA file and groups transcripts and sequences by Gene ID.
    """
    print(f"Reading and parsing FASTA file: {fasta_path}...")
    genes_data = {}
    
    current_protein_id = None
    current_gene_id = None
    current_transcript_id = None
    current_seq_parts = []
    
    count = 0
    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # Save previous record if it exists
                if current_protein_id and current_gene_id and current_transcript_id:
                    seq = "".join(current_seq_parts)
                    if current_gene_id not in genes_data:
                        genes_data[current_gene_id] = []
                    genes_data[current_gene_id].append({
                        "protein_id": current_protein_id,
                        "transcript_id": current_transcript_id,
                        "sequence": seq
                    })
                    count += 1
                    if count % 20000 == 0:
                        print(f"  Parsed {count} protein sequences...")
                
                # Reset and parse new header
                current_seq_parts = []
                parts = line.split()
                # Protein ID (e.g., >ENSP00000481738.1)
                raw_protein = parts[0][1:]
                current_protein_id = raw_protein.split(".")[0]
                
                current_gene_id = None
                current_transcript_id = None
                
                # Find gene and transcript tags in the rest of the header
                for part in parts[1:]:
                    if part.startswith("gene:"):
                        current_gene_id = part.split(":")[1].split(".")[0]
                    elif part.startswith("transcript:"):
                        current_transcript_id = part.split(":")[1].split(".")[0]
            else:
                current_seq_parts.append(line)
        
        # Save the last record
        if current_protein_id and current_gene_id and current_transcript_id:
            seq = "".join(current_seq_parts)
            if current_gene_id not in genes_data:
                genes_data[current_gene_id] = []
            genes_data[current_gene_id].append({
                "protein_id": current_protein_id,
                "transcript_id": current_transcript_id,
                "sequence": seq
            })
            count += 1
            
    print(f"Successfully parsed {count} sequences across {len(genes_data)} unique genes.")
    return genes_data

def build_cache_database(genes_data, db_path):
    """
    Compiles grouped gene transcripts into database entries and writes them to shelve.
    """
    print(f"Creating/updating the shelve cache database at {db_path}...")
    
    # Ensure parent directory exists
    os.makedirs(os.path.dirname(db_path), exist_ok=True)
    
    # Clean up old database files if starting fresh to prevent corruption
    for ext in ['.dat', '.dir', '.bak', '.db']:
        path_to_remove = db_path + ext
        if os.path.exists(path_to_remove):
            os.remove(path_to_remove)
            
    with shelve.open(db_path, flag='c') as db:
        for idx, (gene_id, isoforms) in enumerate(genes_data.items()):
            # Group identical protein sequences to avoid redundant runs
            unique_seqs = {}
            for iso in isoforms:
                seq = iso["sequence"]
                p_id = iso["protein_id"]
                t_id = iso["transcript_id"]
                
                if seq not in unique_seqs:
                    unique_seqs[seq] = {
                        "protein_id": p_id,
                        "transcript_ids": []
                    }
                unique_seqs[seq]["transcript_ids"].append(t_id)
            
            # Construct transcripts_ids list and transcripts_mapping dictionary
            transcripts_ids = []
            transcripts_mapping = {}
            fasta_lines = []
            
            for seq, info in unique_seqs.items():
                p_id = info["protein_id"]
                t_ids = info["transcript_ids"]
                
                transcripts_ids.extend(t_ids)
                
                # Mapping maps representative protein ID to all transcript IDs translating to it
                transcripts_mapping[p_id] = "<br>".join(t_ids)
                
                fasta_lines.append(f">{p_id}\n{seq}")
            
            fasta_content = "\n".join(fasta_lines) + "\n"
            one_isoform = len(unique_seqs) == 1
            
            # Write key-value store entry for the gene
            db[gene_id] = {
                "fasta_content": fasta_content,
                "transcripts_ids": transcripts_ids,
                "transcripts_mapping": transcripts_mapping,
                "one_isoform": one_isoform
            }
            
            if (idx + 1) % 5000 == 0:
                print(f"  Cached data for {idx + 1} genes...")

    print("Success! Ensembl sequence cache is ready.")

def main():
    # Setup paths relative to the script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    root_dir = os.path.dirname(script_dir)
    
    fasta_path = os.path.join(script_dir, "files_genes", "Homo_sapiens.GRCh38.pep.all.fa")
    db_path = os.path.join(root_dir, ".cache_ensembl", "sequences_cache")
    
    # Allow overriding paths via command-line arguments
    if len(sys.argv) > 1:
        fasta_path = sys.argv[1]
    if len(sys.argv) > 2:
        db_path = sys.argv[2]
        
    if not os.path.exists(fasta_path):
        print(f"Error: FASTA file not found at: {fasta_path}")
        print("Please provide the path as an argument: python convert_fasta_to_cache.py <path_to_fasta>")
        sys.exit(1)
        
    genes_data = parse_fasta(fasta_path)
    build_cache_database(genes_data, db_path)

if __name__ == "__main__":
    main()
