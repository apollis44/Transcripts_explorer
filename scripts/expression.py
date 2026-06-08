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

def extracting_tpm_expression_data(ensembl_ids, samples, available_transcripts):
    """
    Given a list of Ensembl transcript IDs, extract their expression data from the TCGA and GTEx datasets.

    Args:
        ensembl_ids (list): A list of Ensembl transcript IDs.
        samples (list): A list of samples to extract expression data for.
        available_transcripts (list): A list of all transcripts available in the dataset.
    Returns:
        pd.DataFrame: A DataFrame containing the expression data for the specified transcripts.
    """

    # Getting the right dataset from Xena
    host = "https://toil.xenahubs.net" # Public hub with both TCGA and GTEx data
    cohort = "TCGA TARGET GTEx" # Cohort name
    dataset = "TcgaTargetGtex_rsem_isoform_tpm" # Dataset name

    available_transcripts_without_version = [t.split('.')[0] for t in available_transcripts]  # Remove version numbers

    # Filter transcripts to those available in the dataset and return their positions
    transcripts_for_expression = [available_transcripts[i] for i in range(len(available_transcripts)) if available_transcripts_without_version[i] in ensembl_ids]
    ensembl_ids_found_in_database = [ensembl_ids[i] for i in range(len(ensembl_ids)) if ensembl_ids[i] in available_transcripts_without_version]
    expression_df = xena.dataset_fetch(host, dataset, samples, sorted(transcripts_for_expression))
    expression_df = pd.DataFrame(expression_df, columns=samples, index=sorted(ensembl_ids_found_in_database))

    return expression_df

def convert_str_to_dict(s):
    s_list = s.split("\t") 
    return {i: s_list[i] for i in range(len(s_list))}

def get_samples_metadata(samples, codes):
    """
    Fetch and process the samples metadata from the TCGA and GTEx datasets.

    Args:
        samples (list): A list of samples to extract metadata for.
        codes (list): A list of codes to convert the metadata to their actual values.
    Returns:
        pd.DataFrame: A DataFrame containing the phenotype metadata.
        pd.DataFrame: A DataFrame containing the survival metadata.
    """

    # Getting the right dataset from Xena
    host = "https://toil.xenahubs.net" # Public hub with both TCGA and GTEx data
    cohort = "TCGA TARGET GTEx" # Cohort name
    metadata_dataset = "TcgaTargetGTEX_phenotype.txt" # Metadata dataset name
    fields = ["_study", "detailed_category"] # Fields to extract

    # Fetch the metadata for the samples
    samples_metadata_phenotype = xena.dataset_fetch(host, metadata_dataset, samples, fields)
    samples_metadata_phenotype = pd.DataFrame(samples_metadata_phenotype, columns=samples, index=fields).T
    samples_metadata_phenotype = samples_metadata_phenotype.map(lambda x: np.nan if x == "NaN" else x) # Create NaN values to remove later
    samples_metadata_phenotype.dropna(inplace=True) # Remove samples with NaN values

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
    common_samples = expression_df.columns.intersection(samples_metadata_phenotype.index)

    # We want to keep only samples from GTEX and TCGA
    metadata_filtered = samples_metadata_phenotype.loc[common_samples]
    metadata_filtered = metadata_filtered[metadata_filtered['_study'].isin(["GTEX", "TCGA"])]

    # We transpose and reindex the expression dataframe
    expression_df = expression_df.loc[:, metadata_filtered.index].T

    # We add the metadata columns to the expression dataframe
    expression_df["tissue_type"] = metadata_filtered["_study"] + " " + metadata_filtered["detailed_category"]
    expression_df["study"] = metadata_filtered["_study"]

    return expression_df

def getting_expression_data(tpm_expression_df, metadata_df, ensembl_ids, mapping):
    common_ids = [x for x in tpm_expression_df.obs_names if x.split('.')[0] in set(ensembl_ids)]
    expression_tpm_df_init = tpm_expression_df[common_ids, :].to_memory().to_df()
    expression_tpm_df_init.index = [x.split('.')[0] for x in expression_tpm_df_init.index]

    if len(expression_tpm_df_init) == 0:
        return None

    # The expression data is: log2(tpm + 0.001), we convert it to tpm
    # Vectorized conversion using NumPy
    vals = expression_tpm_df_init.values.astype(float)
    vals = np.power(2.0, vals) - 0.001
    expression_tpm_df_init = pd.DataFrame(vals, index=expression_tpm_df_init.index, columns=expression_tpm_df_init.columns)

    # We sum the expression values for each isoform
    merged_data = {}
    for isoform in dict.fromkeys(mapping.values()):
        transcripts = isoform.split("<br>")
        common_transcripts = [t for t in transcripts if t in expression_tpm_df_init.index]
        if len(common_transcripts) > 0:
            new_index_name = "<br>".join(common_transcripts)
            merged_data[new_index_name] = expression_tpm_df_init.loc[common_transcripts, :].sum(axis=0)

    if merged_data:
        expression_tpm_df_merged = pd.DataFrame.from_dict(merged_data, orient='index')
    else:
        expression_tpm_df_merged = pd.DataFrame(columns=expression_tpm_df_init.columns)

    # We normalize the expression data to log2(tpm + 1) using vectorized NumPy operation
    vals_merged = expression_tpm_df_merged.values.astype(float)
    expression_normalized_df = pd.DataFrame(np.log2(vals_merged + 1.0), index=expression_tpm_df_merged.index, columns=expression_tpm_df_merged.columns)

    metadata_df_phenotype_init = metadata_df
    expression_normalized_df = merge_expression_and_metadata(expression_normalized_df, metadata_df_phenotype_init)

    return expression_normalized_df

def create_expression_figure_objects(expression_df):

    expression_df = expression_df[expression_df["study"].isin(["GTEX", "TCGA"])]
    dd = expression_df.melt(id_vars=["tissue_type", "study"], var_name="protein", value_name="expression")

    dd = dd.sort_values(by=["study", "tissue_type"])
    proteins = dd.loc[:, "protein"].unique()

    def get_boxplot_stats(group_data):
        """
        Calculates standard boxplot statistics for a given array of data.
        """
        vals = group_data.values
        if len(vals) == 0:
            return pd.Series({
                "median": np.nan,
                "q1": np.nan,
                "q3": np.nan,
                "lowerfence": np.nan,
                "upperfence": np.nan,
                "y": [],
            })

        q1, med, q3 = np.percentile(vals, [25, 50, 75])
        iqr = q3 - q1

        low_limit = q1 - 1.5 * iqr
        high_limit = q3 + 1.5 * iqr

        inside_mask = (vals >= low_limit) & (vals <= high_limit)
        inside_data = vals[inside_mask]

        lowerfence = inside_data.min() if len(inside_data) > 0 else q1
        upperfence = inside_data.max() if len(inside_data) > 0 else q3

        fliers = vals[~inside_mask].tolist()

        return pd.Series({
            "median": med,
            "q1": q1,
            "q3": q3,
            "lowerfence": lowerfence,
            "upperfence": upperfence,
            "y": fliers,
        })

    dd = dd.groupby(['study', 'tissue_type', 'protein'], sort=False)['expression'].apply(get_boxplot_stats).unstack().reset_index()

    return dd