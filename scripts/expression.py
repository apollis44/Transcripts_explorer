import xenaPython as xena
import numpy as np
import pandas as pd
import ssl
import seaborn as sns
import matplotlib.pyplot as plt
import os
import json

# Bypass SSL certificate verification for Xena Hub
ssl._create_default_https_context = ssl._create_unverified_context

def extracting_count_expression_data(ensembl_ids):
    """
    Given a list of Ensembl transcript IDs, extract their expression data from the TCGA and GTEx datasets.

    Args:
        ensembl_ids (list): A list of Ensembl transcript IDs.
    Returns:
        pd.DataFrame: A DataFrame containing the expression data for the specified transcripts.
    """

    # Getting the right dataset from Xena
    host = "https://toil.xenahubs.net" # Public hub with both TCGA and GTEx data
    cohort = "TCGA TARGET GTEx" # Cohort name
    samples = xena.cohort_samples(host, cohort, None)
    dataset = "TcgaTargetGtex_expected_count" # Dataset name

    # We extract the list of all transcripts ids available in the dataset
    available_transcripts = xena.dataset_field(host, dataset)
    available_transcripts_without_version = [t.split('.')[0] for t in available_transcripts]  # Remove version numbers

    # Filter transcripts to those available in the dataset and return their positions
    transcripts_for_expression = [available_transcripts[i] for i in range(len(available_transcripts)) if available_transcripts_without_version[i] in ensembl_ids]
    ensembl_ids_found_in_database = [ensembl_ids[i] for i in range(len(ensembl_ids)) if ensembl_ids[i] in available_transcripts_without_version]
    for transcript in ensembl_ids:
        if transcript not in ensembl_ids_found_in_database:
            print(f"\033[31mWarning: Transcript {transcript} not found in the database.\033[0m")
    expression_df = xena.dataset_fetch(host, dataset, samples, sorted(transcripts_for_expression))
    expression_df = pd.DataFrame(expression_df, columns=samples, index=sorted(ensembl_ids_found_in_database))

    return expression_df

def extracting_tpm_expression_data(ensembl_ids):
    """
    Given a list of Ensembl transcript IDs, extract their expression data from the TCGA and GTEx datasets.

    Args:
        ensembl_ids (list): A list of Ensembl transcript IDs.
    Returns:
        pd.DataFrame: A DataFrame containing the expression data for the specified transcripts.
    """

    # Getting the right dataset from Xena
    host = "https://toil.xenahubs.net" # Public hub with both TCGA and GTEx data
    cohort = "TCGA TARGET GTEx" # Cohort name
    samples = xena.cohort_samples(host, cohort, None)
    dataset = "TcgaTargetGtex_rsem_isoform_tpm" # Dataset name

    # We extract the list of all transcripts ids available in the dataset
    available_transcripts = xena.dataset_field(host, dataset)
    available_transcripts_without_version = [t.split('.')[0] for t in available_transcripts]  # Remove version numbers

    # Filter transcripts to those available in the dataset and return their positions
    transcripts_for_expression = [available_transcripts[i] for i in range(len(available_transcripts)) if available_transcripts_without_version[i] in ensembl_ids]
    ensembl_ids_found_in_database = [ensembl_ids[i] for i in range(len(ensembl_ids)) if ensembl_ids[i] in available_transcripts_without_version]
    for transcript in ensembl_ids:
        if transcript not in ensembl_ids_found_in_database:
            print(f"\033[31mWarning: Transcript {transcript} not found in the database.\033[0m")
    expression_df = xena.dataset_fetch(host, dataset, samples, sorted(transcripts_for_expression))
    expression_df = pd.DataFrame(expression_df, columns=samples, index=sorted(ensembl_ids_found_in_database))

    return expression_df

def convert_str_to_dict(s):
    s_list = s.split("\t") 
    return {i: s_list[i] for i in range(len(s_list))}

def get_samples_metadata():
    """
    Fetch and process the samples metadata from the TCGA and GTEx datasets.

    Returns:
        pd.DataFrame: A DataFrame containing the phenotype metadata.
        pd.DataFrame: A DataFrame containing the survival metadata.
    """

    # Getting the right dataset from Xena
    host = "https://kidsfirst.xenahubs.net" # Public hub with both TCGA and GTEx data
    cohort = "TCGA TARGET GTEx KidsFirst" # Cohort name
    samples = xena.cohort_samples(host, cohort, None)
    metadata_dataset = "TCGA_target_GTEX_KF/phenotype.txt" # Metadata dataset name
    fields = ["_study", "main_category"] # Fields to extract

    # Fetch the metadata for the samples
    samples_metadata_phenotype = xena.dataset_fetch(host, metadata_dataset, samples, fields)
    samples_metadata_phenotype = pd.DataFrame(samples_metadata_phenotype, columns=samples, index=fields).T
    samples_metadata_phenotype = samples_metadata_phenotype.map(lambda x: np.nan if x == "NaN" else x) # Create NaN values to remove later
    samples_metadata_phenotype.dropna(inplace=True) # Remove samples with NaN values

    # The dataset contains codes for categorical fields, we convert them to their actual values
    # To do so, we fetch the codes from Xena and create a mapping dictionary
    codes = xena.field_codes(host, metadata_dataset, fields)

    for i, field in enumerate(fields):
        codes_field = codes[i]['code']
        codes_field_dict = convert_str_to_dict(codes_field)
        samples_metadata_phenotype[field] = samples_metadata_phenotype[field].astype(object)
        for i in range(len(samples_metadata_phenotype)):
            samples_metadata_phenotype.iloc[i, samples_metadata_phenotype.columns.get_loc(field)] = codes_field_dict[int(samples_metadata_phenotype.iloc[i][field])]

    return samples_metadata_phenotype

def merge_expression_and_metadata(expression_df, samples_metadata_phenotype):
    """
    Merge the expression data with the samples metadata.

    Args:
        expression_df (pd.DataFrame): DataFrame containing expression data.
        samples_metadata_phenotype (pd.DataFrame): DataFrame containing phenotype metadata.

    Returns:
        pd.DataFrame: Merged DataFrame containing both expression data and metadata.
    """

    # We keep only samples that are in the datasets
    common_samples = expression_df.columns.intersection(samples_metadata_phenotype.index).tolist()
    expression_df = expression_df.loc[:,common_samples].T
    samples_metadata_phenotype = samples_metadata_phenotype.loc[common_samples]

    # We want to keep only samples from GTEX and TCGA
    samples_to_keep = samples_metadata_phenotype.loc[samples_metadata_phenotype['_study'].isin(["GTEX", "TCGA"]),:].index.tolist()
    expression_df = expression_df.loc[samples_to_keep]
    samples_metadata_phenotype = samples_metadata_phenotype.loc[samples_to_keep]

    # We add the metadata columns to the expression dataframe
    expression_df["cancer_type"] = samples_metadata_phenotype.loc[expression_df.index, "main_category"]
    expression_df["study"] = samples_metadata_phenotype.loc[expression_df.index, "_study"]

    return expression_df

def extract_normalization_values(expression_count_df_init, expression_tpm_df_init, transcripts_length):
    """
    Extract the normalization values from the expression data.

    Args:
        expression_count_df_init (pd.DataFrame): DataFrame containing the count expression data.
        expression_tpm_df_init (pd.DataFrame): DataFrame containing the TPM expression data.
        transcripts_length (dict): Dictionary containing the length of each transcript.

    Returns:
        np.ndarray: Array containing the normalization values.
    """
    # The values in the dataset are log2(counts + 1) and log2(tpm + 0.001)
    expression_count_numpy_init = expression_count_df_init.to_numpy(dtype=np.float64)
    expression_tpm_numpy_init = expression_tpm_df_init.to_numpy(dtype=np.float64)

    expression_count_numpy_init = np.power(2, expression_count_numpy_init) - 1
    expression_tpm_numpy_init = np.power(2, expression_tpm_numpy_init) - 0.001

    # Replace negative TPM values with 0
    expression_tpm_numpy_init[expression_tpm_numpy_init < 0] = 0

    # TPM is counts / length / normalization_value
    transcripts_length_matched = [transcripts_length[transcript] for transcript in expression_count_df_init.index]
    lengths_vector = np.array(transcripts_length_matched)[:, np.newaxis]
    expression_tpm_numpy_init = expression_tpm_numpy_init * lengths_vector

    normalization_values = np.divide(
        expression_count_numpy_init, 
        expression_tpm_numpy_init, 
        out=np.zeros_like(expression_count_numpy_init), 
        where=expression_tpm_numpy_init != 0
    )
    masked_data = np.ma.masked_equal(normalization_values, 0)
    normalization_values = masked_data.mean(axis=0).filled(0)

    return normalization_values

def normalize_expression(expression_count_df, normalization_values, transcripts_length):
    """
    Normalize the expression data.

    Args:
        expression_count_df (pd.DataFrame): DataFrame containing the count expression data.
        normalization_values (np.ndarray): Array containing the normalization values.
        transcripts_length (dict): Dictionary containing the length of each transcript.

    Returns:
        pd.DataFrame: DataFrame containing the normalized expression data.
    """
    expression_count_numpy = expression_count_df.to_numpy(dtype=np.float64)
    expression_count_numpy = np.divide(
        expression_count_numpy, 
        normalization_values, 
        out=np.zeros_like(expression_count_numpy), # Fill output with 0s initially
        where=normalization_values != 0             # Only perform division where divisor is not 0
    )
    
    transcripts_length_matched = [sum(transcripts_length[t_id] for t_id in isoform.split("<br>")) for isoform in expression_count_df.index]

    lengths_vector = np.array(transcripts_length_matched, dtype=np.float64)[:, np.newaxis]

    expression_count_numpy = expression_count_numpy / lengths_vector

    # We go back to log2(counts + 1)
    expression_count_numpy = np.log2(expression_count_numpy + 1)

    expression_count_df = pd.DataFrame(expression_count_numpy, columns=expression_count_df.columns, index=expression_count_df.index)

    return expression_count_df 

def getting_expression_data(ensembl_ids, transcripts_length, mapping):
    expression_count_df_init = extracting_count_expression_data(ensembl_ids)
    expression_tpm_df_init = extracting_tpm_expression_data(ensembl_ids)

    # Match columns and rows
    expression_count_df_init = expression_count_df_init.reindex_like(expression_tpm_df_init)

    normalization_values = extract_normalization_values(expression_count_df_init, expression_tpm_df_init, transcripts_length)

    # The expression data is: log2(counts + 1), we convert it to counts
    # We first convert values to float to avoid errors
    expression_count_df_init = expression_count_df_init.map(lambda x: float(x))
    expression_count_df_init = expression_count_df_init.map(lambda x: np.power(2, x) - 1)

    expression_count_df_merged = pd.DataFrame(columns=expression_count_df_init.columns)
    # We sum the expression values for each isoform
    for isoform in dict.fromkeys(mapping.values()):
        transcripts = isoform.split("<br>")
        common_transcripts = list(set(transcripts) & set(expression_count_df_init.index))
        if len(common_transcripts) > 0:
            new_index_name = "<br>".join(common_transcripts)
            expression_count_df_merged.loc[new_index_name,:] = expression_count_df_init.loc[common_transcripts,:].sum(axis=0)

    # We normalize the expression data
    expression_normalized_df = normalize_expression(expression_count_df_merged, normalization_values, transcripts_length)

    metadata_df_phenotype_init = get_samples_metadata()
    expression_normalized_df = merge_expression_and_metadata(expression_normalized_df, metadata_df_phenotype_init)

    return expression_normalized_df

def create_expression_figure_objects(expression_df):
    print(f"Creating expression figure objects...")

    expression_df = expression_df[expression_df["study"].isin(["GTEX", "TCGA"])]
    dd = expression_df.melt(id_vars=["cancer_type", "study"], var_name="protein", value_name="expression")

    dd = dd.sort_values(by=["study", "cancer_type"])
    proteins = dd.loc[:, "protein"].unique()

    def get_boxplot_stats(group_data):
        """
        Calculates standard boxplot statistics for a given array of data.
        """
        # 1. Basic Stats
        med = group_data.median()
        q1 = group_data.quantile(0.25)
        q3 = group_data.quantile(0.75)
        iqr = q3 - q1
        
        # 2. Whiskers (1.5 * IQR rule)
        low_limit = q1 - 1.5 * iqr
        high_limit = q3 + 1.5 * iqr
        
        # Identify data points within the limits
        inside_data = group_data[(group_data >= low_limit) & (group_data <= high_limit)]
        
        # The whisker ends are the actual min/max of the data INSIDE the limits
        # (Matplotlib convention: whiskers don't extend to empty space)
        lowerfence = inside_data.min() if not inside_data.empty else q1
        upperfence = inside_data.max() if not inside_data.empty else q3
        
        # 3. Outliers (Fliers)
        fliers = group_data[(group_data < lowerfence) | (group_data > upperfence)].tolist()
        
        return pd.Series({
            "median": med,
            "q1": q1,
            "q3": q3,
            "lowerfence": lowerfence,
            "upperfence": upperfence,
            "y": fliers,
        })
    
    dd = dd.groupby(['study', 'cancer_type', 'protein'], sort=False)['expression'].apply(get_boxplot_stats).unstack().reset_index()

    return dd