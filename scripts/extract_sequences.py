from matplotlib.pylab import result_type
import requests
import json
import time

def fetch_transcripts(ensembl_id, session):
    """
    Fetches the transcripts for a given Ensembl gene ID from the Ensembl REST API.

    Args:
        ensembl_id (str): The gene Ensembl ID (e.g., ENSG00000156738).

    Returns:
        str: The JSON response containing isoform information, or None if an error occurs.
    """
    server = "https://rest.ensembl.org"

    # The endpoint for retrieving isoforms by gene ID.
    endpoint = f"/lookup/id/{ensembl_id}"
    
    # Headers to specify that we want the response in JSON format.
    headers = {"Content-Type": "application/json"}
    
    # Parameters to retrieve expanded information including transcripts.
    params = {"expand": "1"}

    delay = 2
    for attempt in range(10):
        try:
            # Make the GET request to the Ensembl API.
            response = session.get(server + endpoint, headers=headers, params=params, timeout=45)
            
            if response.status_code == 429:
                retry_after = response.headers.get("Retry-After")
                wait_time = int(retry_after) if retry_after and retry_after.isdigit() else delay
                print(f"Rate limited (429). Waiting {wait_time} seconds before retry...")
                time.sleep(wait_time)
                delay = min(delay * 2, 60)
                continue

            # Raise an exception for bad status codes (4xx or 5xx).
            response.raise_for_status()

            # If the request was successful, return the json output.
            res_json = response.json()

            transcripts = []
            for transcript in res_json['Transcript']:
                if transcript["biotype"] == "protein_coding":
                    transcripts.append(transcript)
            return transcripts

        except Exception as e:
            print(f"Error fetching transcripts for {ensembl_id} (attempt {attempt + 1}/10): {e}")
            if attempt == 9:
                raise RuntimeError(f"Failed to fetch transcripts for {ensembl_id} after 10 attempts: {e}") from e
            time.sleep(delay)
            delay = min(delay * 2, 60)

def fetch_protein_sequence(ensembl_ids, session):
    """
    Fetches the protein sequence for a given Ensembl ID from the Ensembl REST API.

    Args:
        ensembl_ids (list): list of ensembl ids of the proteins to fetch the sequence from. 

    Returns:
        list: list of protein sequences in FASTA format, or None if an error occurs.
    """
    server = "https://rest.ensembl.org"
    
    # The endpoint for retrieving a sequence by its ID.
    # We specify the 'protein' type to ensure we get the protein sequence.
    endpoint = "/sequence/id/"

    data = {"ids": ensembl_ids}
    
    # Headers to specify that we want the response in FASTA format.
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    
    # Parameters to specify the type of sequence to retrieve.
    params = {"type": "protein"}

    delay = 2
    for attempt in range(10):
        try:
            # Make the POST request to the Ensembl API.
            response = session.post(server+endpoint, headers=headers, data=json.dumps(data), params=params, timeout=45)
            
            if response.status_code == 429:
                retry_after = response.headers.get("Retry-After")
                wait_time = int(retry_after) if retry_after and retry_after.isdigit() else delay
                print(f"Rate limited (429). Waiting {wait_time} seconds before retry...")
                time.sleep(wait_time)
                delay = min(delay * 2, 60)
                continue

            # Raise an exception for bad status codes (4xx or 5xx).
            response.raise_for_status()

            res_json = response.json()

            protein_sequences = {item["query"]: item["seq"] for item in res_json}
            
            # If the request was successful, return the sequence text.
            return protein_sequences

        except Exception as e:
            print(f"Error fetching protein sequences (attempt {attempt + 1}/10): {e}")
            if attempt == 9:
                raise RuntimeError(f"Failed to fetch protein sequences after 10 attempts: {e}") from e
            time.sleep(delay)
            delay = min(delay * 2, 60)

def fetch_cdna_length(ensembl_ids, session):
    """
    Fetches the cDNA (coding sequence) for a given Ensembl ID from the Ensembl REST API.

    Args:
        ensembl_ids (list): list of ensembl ids of the cdnas to fetch the sequence from. 
    Returns:
        dict: dictionary of cdna sequences in FASTA format, or None if an error occurs.
    """
    server = "https://rest.ensembl.org"
    
    # The endpoint for retrieving a sequence by its ID.
    # We specify the 'transcript' type to ensure we get the transcript sequence.
    endpoint = "/sequence/id/"

    data = {"ids": ensembl_ids}
    
    # Headers to specify that we want the response in FASTA format.
    headers = {"Content-Type": "application/json", "Accept": "application/json"}

    # Parameters to specify the type of sequence to retrieve.
    params = {"type": "cdna"}
    delay = 2
    for attempt in range(10):
        try:
            # Make the POST request to the Ensembl API.
            response = session.post(server + endpoint, headers=headers, data=json.dumps(data), params=params, timeout=45)
            
            if response.status_code == 429:
                retry_after = response.headers.get("Retry-After")
                wait_time = int(retry_after) if retry_after and retry_after.isdigit() else delay
                print(f"Rate limited (429). Waiting {wait_time} seconds before retry...")
                time.sleep(wait_time)
                delay = min(delay * 2, 60)
                continue

            # Raise an exception for bad status codes (4xx or 5xx).
            response.raise_for_status()

            res_json = response.json()
            
            # If the request was successful, return the sequence text.
            return {item["query"]: len(item["seq"]) for item in res_json}

        except Exception as e:
            print(f"Error fetching cDNA length (attempt {attempt + 1}/10): {e}")
            if attempt == 9:
                raise RuntimeError(f"Failed to fetch cDNA length after 10 attempts: {e}") from e
            time.sleep(delay)
            delay = min(delay * 2, 60)
