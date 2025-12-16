from matplotlib.pylab import result_type
import requests
import time
import sys
import json

def fetch_transcripts(ensembl_id):
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

    print(f"Fetching transcripts for gene with Ensembl ID: {ensembl_id}...")

    try:
        # Make the GET request to the Ensembl API.
        response = requests.get(server + endpoint, headers=headers, params=params)
        
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

def fetch_protein_sequence(ensembl_ids):
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

    print(f"Fetching protein sequence for Ensembl IDs: {' '.join(ensembl_ids)}")

    try:
        # Make the GET request to the Ensembl API.
        response = requests.post(server+endpoint, headers=headers, data=json.dumps(data), params=params)
        
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

def fetch_cds(ensembl_id):
    """
    Fetches the CDS (coding sequence) for a given Ensembl ID from the Ensembl REST API.

    Args:
        ensembl_id (str): The Ensembl transcript or gene ID.
                          (e.g., ENST00000380152).
    Returns:
        str: The transcript sequence in FASTA format, or None if an error occurs.
    """
    server = "https://rest.ensembl.org"
    
    # The endpoint for retrieving a sequence by its ID.
    # We specify the 'transcript' type to ensure we get the transcript sequence.
    endpoint = f"/sequence/id/{ensembl_id}"
    
    # Headers to specify that we want the response in FASTA format.
    headers = {"Content-Type": "text/plain"}

    # Parameters to specify the type of sequence to retrieve.
    params = {"type": "cds"}

    print(f"Fetching CDS for Ensembl ID: {ensembl_id}...")

    try:
        # Make the GET request to the Ensembl API.
        response = requests.get(server + endpoint, headers=headers, params=params)
        
        # Raise an exception for bad status codes (4xx or 5xx).
        response.raise_for_status()
        
        # If the request was successful, return the sequence text.
        return response.text

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
    
def fetch_sequence_with_introns_information(ensembl_id):
    """
    Fetches the genomic sequence for a given Ensembl ID from the Ensembl REST API.

    Args:
        ensembl_id (str): The Ensembl transcript or gene ID.

    Returns:
        str: The transcript sequence in FASTA format, or None if an error occurs. 
             The sequence is in capital letters for exons and lowercase for introns.
    """
    server = "https://rest.ensembl.org"
    
    # The endpoint for retrieving a sequence by its ID.
    # We specify the 'transcript' type to ensure we get the transcript sequence.
    endpoint = f"/sequence/id/{ensembl_id}"
    
    # Headers to specify that we want the response in FASTA format.
    headers = {"Content-Type": "text/plain"}

    # Parameters to specify the type of sequence to retrieve.
    params = {"mask_feature": "1"}

    print(f"Fetching sequence with introns information for Ensembl ID: {ensembl_id}...")

    try:
        # Make the GET request to the Ensembl API.
        response = requests.get(server + endpoint, headers=headers, params=params)
        
        # Raise an exception for bad status codes (4xx or 5xx).
        response.raise_for_status()
        
        # If the request was successful, return the sequence text.
        return response.text

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

def fetch_canonical_uniprot_id(gene_name):
    """
    Fetches the canonical UniProt ID for a given gene name from the UniProt REST API.

    Args:
        gene_name (str): The gene name.

    Returns:
        str: The canonical UniProt ID, or None if an error occurs.
    """
    uniprot_api = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f"gene:{gene_name} AND organism_id:9606", # 9606 = Human
        "fields": "accession,id",
        "format": "json"
    }
    
    try:
        # Make the GET request to the UniProt API.
        response = requests.get(uniprot_api, params=params)
        
        # Raise an exception for bad status codes (4xx or 5xx).
        response.raise_for_status()
        
        # If the request was successful, return the canonical UniProt ID.
        print("Canonical UniProt ID: ", response.json()["results"][0]["primaryAccession"])
        return response.json()["results"][0]["primaryAccession"]
    except requests.exceptions.HTTPError as http_err:
        # Handle HTTP errors (e.g., 404 Not Found, 500 Server Error).
        print(f"\033[31mHTTP error occurred: {http_err}\033[0m")
        print(f"Status Code: {response.status_code}")
        print("Please check if the canonical Ensembl ID is correct.")
        # The response content might contain a more specific error message from the server.
        print(f"Server response: {response.content.decode()}")
        return None
    except requests.exceptions.RequestException as req_err:
        # Handle other request errors (e.g., network issues).
        print(f"\033[31mAn error occurred: {req_err}\033[0m")
        return None
    
def align_sequences(sequences, email):
    """
    Aligns the given sequences using the Clustal Omega service from EBI.

    Args:
        sequences (str): The sequences in FASTA format to be aligned.
        email (str): A valid email address for identification.

    Returns:
        str: The aligned sequences in FASTA format, or None if an error occurs.
    """
    # --- 1. DEFINE CORRECT BASE URL AND PARAMETERS ---

    # The correct base URL for the EBI REST API
    BASE_URL = "https://www.ebi.ac.uk/Tools/services/rest/clustalo"
    RUN_URL = f"{BASE_URL}/run"
    STATUS_URL = f"{BASE_URL}/status"
    RESULT_URL = f"{BASE_URL}/result"

    # --- 2. FUNCTION TO SUBMIT JOB ---
    def submit_job(sequences, email):
        """
        Submits sequences to Clustal Omega and returns the job ID.

        Args:
            sequences (str): The sequences in FASTA format.
            email (str): A valid email address for identification.

        Returns:
            str: The job ID if submission is successful, None otherwise.
        """
        print("Submitting alignment job...")
        # Parameters for the POST request. Keys must be exact.
        params = {
            'email': email,
            'sequence': sequences,
        }

        try:
            # Use POST for submission
            response = requests.post(RUN_URL, data=params)
            response.raise_for_status()  # Raises an exception for bad status codes (4xx or 5xx)

            # The correct API returns the job ID as plain text
            job_id = response.text
            print(f"Successfully submitted. Job ID: {job_id}")
            return job_id
        except requests.exceptions.RequestException as e:
            print(f"\033[31mError submitting job: {e}\033[0m")
            # Print server response if available
            if e.response:
                print(e.response.text)
            return None

    # --- 3. FUNCTION TO CHECK STATUS ---
    def check_status(job_id):
        """
        Checks the status of a running job.

        Args:
            job_id (str): The job ID to check.

        Returns:
            str: The status of the job.
        """
        try:
            status_endpoint = f"{STATUS_URL}/{job_id}"
            response = requests.get(status_endpoint)
            response.raise_for_status()
            return response.text
        except requests.exceptions.RequestException as e:
            print(f"\033[31mError checking status for job {job_id}: {e}\033[0m")
            return "ERROR"

    # --- 4. FUNCTION TO GET RESULT ---
    def get_result(job_id, result_type='fa'):
        """
        Retrieves the result for a finished job.
        Args:
            job_id (str): The job ID.
            result_type (str): The type of result to fetch (default is 'fa' for FASTA). 

        Returns:
            str: The result content.
        """
        print("Fetching alignment result...")
        try:
            result_endpoint = f"{RESULT_URL}/{job_id}/{result_type}"
            response = requests.get(result_endpoint)
            response.raise_for_status()
            return response.text
        except requests.exceptions.RequestException as e:
            print(f"\033[31mError fetching result for job {job_id}: {e}\033[0m")
            sys.exit(1)

    # Step 1: Submit the job
    job_id = submit_job(sequences, email)

    # Step 2: Poll for status until finished
    status = ""
    while status != "FINISHED":
        status = check_status(job_id)
        print(f"Job status: {status}")
        if status == "FINISHED":
            break
        elif status in ["ERROR", "NOT_FOUND", "FAILURE"]:
            print("\033[31mJob failed or could not be found.\033[0m")
            sys.exit(1)
        # Wait for 10 seconds before checking again
        time.sleep(10)

    # Step 3: Get and print the result
    alignment = get_result(job_id)

    return alignment