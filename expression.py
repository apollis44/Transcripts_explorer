import xenaPython as xena
import numpy as np
import pandas as pd
import ssl
import os
import matplotlib.pyplot as plt
import seaborn as sns

# Bypass SSL certificate verification for Xena Hub
ssl._create_default_https_context = ssl._create_unverified_context

def extracting_expression_data(ensembl_ids):
    """
    Given a list of Ensembl transcript IDs, extract their expression data from the TCGA and GTEx datasets.
    """
    # Getting the right dataset from Xena
    host = "https://kidsfirst.xenahubs.net" # Public hub with both TCGA and GTEx data
    cohort = "TCGA TARGET GTEx KidsFirst" # Cohort name
    samples = xena.cohort_samples(host, cohort, None)
    dataset = "TCGA_target_GTEX_KF/rsem.isoforms_TPM.txt" # Dataset name

    # We extract the list of all transcripts ids available in the dataset
    available_transcripts = xena.dataset_field(host, dataset)
    available_transcripts_without_version = [t.split('.')[0] for t in available_transcripts]  # Remove version numbers

    # Filter transcripts to those available in the dataset and return their positions
    transcripts_for_expression = [available_transcripts[i] for i in range(len(available_transcripts)) if available_transcripts_without_version[i] in ensembl_ids]
    ensembl_ids_found_in_database = [ensembl_ids[i] for i in range(len(ensembl_ids)) if ensembl_ids[i] in available_transcripts_without_version]
    
    for transcript in ensembl_ids:
        if transcript not in ensembl_ids_found_in_database:
            print(f"Warning: Transcript {transcript} not found in the database.")
            
    expression_df = xena.dataset_fetch(host, dataset, samples, sorted(transcripts_for_expression))
    expression_df = pd.DataFrame(expression_df, columns=samples, index=sorted(ensembl_ids_found_in_database))

    return expression_df

def convert_str_to_dict(s):
    s_list = s.split("\t") 
    return {i: s_list[i] for i in range(len(s_list))}

def get_samples_metadata():
    """
    Fetch and process the samples metadata from the TCGA and GTEx datasets.
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
    codes = xena.field_codes(host, metadata_dataset, fields)

    for i, field in enumerate(fields):
        codes_field = codes[i]['code']
        codes_field_dict = convert_str_to_dict(codes_field)
        samples_metadata_phenotype[field] = samples_metadata_phenotype[field].astype(object)
        # Avoid row-wise loop for performance, though leaving consistent with original logic for safety
        for j in range(len(samples_metadata_phenotype)):
             # Using efficient iloc access
             val = samples_metadata_phenotype.iloc[j, samples_metadata_phenotype.columns.get_loc(field)]
             samples_metadata_phenotype.iloc[j, samples_metadata_phenotype.columns.get_loc(field)] = codes_field_dict[int(val)]

    return samples_metadata_phenotype

def merge_expression_and_metadata(expression_df, samples_metadata_phenotype):
    """
    Merge the expression data with the samples metadata.
    """
    # We keep only samples that are in the three datasets
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

def ensure_dir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def get_expression_data(ensembl_ids, output_dir):
    ensure_dir(output_dir)
    file_path = f"{output_dir}/TCGA_GTEx_expression_data.csv"
    
    if os.path.exists(file_path):
        print("Loading expression data from file...")
        return pd.read_csv(file_path, index_col=0)
    
    print("Fetching expression data from Xena...")
    expression_df = extracting_expression_data(ensembl_ids)
    samples_metadata_phenotype = get_samples_metadata()
    expression_df = merge_expression_and_metadata(expression_df, samples_metadata_phenotype)
    
    expression_df.to_csv(file_path)
    print(f"Saved TCGA_GTEx_expression_data.csv to {output_dir}")
    return expression_df

def plot_expression(expression_df, mapping, type_data="TCGA", specific_cancer_type=None):
    print(f"Plotting {type_data} expression...")

    if specific_cancer_type:
        expression_df = expression_df[expression_df["cancer_type"] == specific_cancer_type]
        number_transcripts = len(expression_df.columns) - 2
        length_figure = number_transcripts + 1
    else:
        length_figure = 15
        
    expression_df = expression_df[expression_df["study"] == type_data]
    nb_cancer_types = len(expression_df["cancer_type"].unique())
    
    if nb_cancer_types == 0:
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "No data available", ha='center', va='center')
        return fig

    # Calculate layout
    cols = 2 if nb_cancer_types > 1 else 1
    rows = nb_cancer_types // 2 + 1 if nb_cancer_types > 1 else 1
    
    fig, axs = plt.subplots(rows, cols, figsize=(length_figure, 5*rows), sharey=True)
    
    dd = expression_df.melt(id_vars=["cancer_type", "study"], var_name="transcript", value_name="expression")

    cancer_types = expression_df["cancer_type"].unique()
    
    # Handle single subplot case
    if nb_cancer_types == 1:
        axs = np.array([axs])

    axs_flat = axs.flatten()

    for i, cancer_type in enumerate(cancer_types):
        current_plot = axs_flat[i]
        dd_cancer = dd[dd["cancer_type"] == cancer_type]
        dd_cancer = dd_cancer.assign(
            protein=dd_cancer["transcript"].map(
                lambda x: next((mapping[name] for name in mapping.keys() if x in name), 
                            "Unknown")
            )
        )
        dd_cancer = dd_cancer.sort_values(by="protein")

        sns.boxplot(data=dd_cancer, 
                    x="transcript",
                    y="expression", 
                    hue="protein", ax=current_plot)
        
        if current_plot.legend_:
            current_plot.legend_.remove()
            
        current_plot.set_ylabel("log2(TPM)")
        current_plot.set_xlabel("")
        current_plot.tick_params(axis='x', rotation=90)
        current_plot.set_title(cancer_type)

    # Hide unused subplots
    for j in range(i + 1, len(axs_flat)):
        axs_flat[j].axis('off')

    # Add legend
    if nb_cancer_types > 1:
        handles, labels = axs_flat[0].get_legend_handles_labels()
        if handles:
            fig.legend(handles, labels, bbox_to_anchor=(1.01, 1), loc='upper left', title="Protein Isoform")
            
    plt.tight_layout()
    return fig
