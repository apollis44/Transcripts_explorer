# %%

import gc
import fastparquet
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix

# --- Configuration ---
parquet_path = "./scripts/files_genes/TcgaTargetGtex_rsem_isoform_tpm.parquet"
output_h5ad_path = "./scripts/files_genes/transcript_counts_sparse.h5ad"

# 1. Open the Parquet file to inspect metadata
print("Reading file metadata via Fastparquet...")
pf = fastparquet.ParquetFile(parquet_path)

# Extract column names (Fastparquet stores them in pf.columns)
all_columns = pf.columns
transcript_id_col = all_columns[0]  # Assuming first column is the ID/Name
sample_names = all_columns[1:]       # Rest are sample columns

print(f"Found {len(sample_names)} samples.")

transcript_ids = []

# 2. Iterate through the Parquet file by Row Groups
print("Streaming file row groups into sparse matrix chunks...")
for df_chunk in pf.iter_row_groups():
    
    # Fastparquet yields a standard pandas DataFrame for each row group
    
    # Extract transcript IDs for the current chunk
    transcript_ids_chunk = df_chunk.index.tolist()
    transcript_ids_chunk = [transcript_id.split(".")[0] for transcript_id in transcript_ids_chunk]  # Remove version suffix if present
    transcript_ids.extend(df_chunk.index.tolist())
    
    # Isolate count data
    counts_only = df_chunk.values
    
    # Convert directly into a Compressed Sparse Row (CSR) matrix (using Float32 for memory efficiency)
    sparse_chunk = csr_matrix(counts_only, dtype="float32")
    
    # Clean up loop variables to protect RAM
    del df_chunk, counts_only
    gc.collect()

# 4. Construct the AnnData Object
print("Building AnnData object...")
adata = ad.AnnData(
    X=sparse_chunk,
    obs=pd.DataFrame(index=transcript_ids),       # Rows = Samples
    var=pd.DataFrame(index=sample_names)      # Columns = Transcripts
)

# 5. Write out the compressed h5ad file
print(f"Saving to {output_h5ad_path}...")
adata.write_h5ad(output_h5ad_path, compression="gzip")

del adata
gc.collect()
print("Done! Sparse AnnData built successfully using Fastparquet backend.")