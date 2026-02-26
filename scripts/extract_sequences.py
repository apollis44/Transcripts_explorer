from matplotlib.pylab import result_type
import requests
import json

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

    try:
        # Make the GET request to the Ensembl API.
        response = session.get(server + endpoint, headers=headers, params=params)
        
        # Raise an exception for bad status codes (4xx or 5xx).
        response.raise_for_status()

        # If the request was successful, return the json output.
        response = response.json()

        transcripts = []
        for transcript in response['Transcript']:
            if transcript["biotype"] == "protein_coding":
                transcripts.append(transcript)
        return transcripts

    except requests.exceptions.HTTPError as http_err:
        # Handle HTTP errors (e.g., 404 Not Found, 500 Server Error).
        print(f"\033[31mHTTP error occurred: {http_err}\033[0m")
        print(f"Status Code: {response.status_code}")
        print("Please check if the Ensembl ID is correct and corresponds to a gene.")
        # The response content might contain a more specific error message from the server.
        print(f"Server response: {response.content.decode()}")
        return None
    except requests.exceptions.RequestException as req_err:
        # Handle other request errors (e.g., network issues).
        print(f"\033[31mAn error occurred: {req_err}\033[0m")
        return None

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

    try:
        # Make the GET request to the Ensembl API.
        response = session.post(server+endpoint, headers=headers, data=json.dumps(data), params=params)
        
        # Raise an exception for bad status codes (4xx or 5xx).
        response.raise_for_status()

        response = response.json()

        protein_sequences = {item["query"]: item["seq"] for item in response}
        
        # If the request was successful, return the sequence text.
        return protein_sequences

    except requests.exceptions.HTTPError as http_err:
        # Handle HTTP errors (e.g., 404 Not Found, 500 Server Error).
        print(f"\033[31mHTTP error occurred: {http_err}\033[0m")
        print(f"Status Code: {response.status_code}")
        print("Please check if the Ensembl ID is correct and corresponds to a protein.")
        # The response content might contain a more specific error message from the server.
        print(f"Server response: {response.content.decode()}")
        return None
    except requests.exceptions.RequestException as req_err:
        # Handle other request errors (e.g., network issues).
        print(f"\033[31mAn error occurred: {req_err}\033[0m")
        return None

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
    try:
        # Make the GET request to the Ensembl API.
        response = session.post(server + endpoint, headers=headers, data=json.dumps(data), params=params)
        
        # Raise an exception for bad status codes (4xx or 5xx).
        response.raise_for_status()

        response = response.json()
        
        # If the request was successful, return the sequence text.
        return {item["query"]: len(item["seq"]) for item in response}

    except requests.exceptions.HTTPError as http_err:
        # Handle HTTP errors (e.g., 404 Not Found, 500 Server Error).
        print(f"\033[31mHTTP error occurred: {http_err}\033[0m")
        print(f"Status Code: {response.status_code}")
        print("Please check if the Ensembl ID is correct and corresponds to a transcript.")
        # The response content might contain a more specific error message from the server.
        print(f"Server response: {response.content.decode()}")
        return None
    except requests.exceptions.RequestException as req_err:
        # Handle other request errors (e.g., network issues).
        print(f"\033[31mAn error occurred: {req_err}\033[0m")
        return None
