# ---
# jupyter:
#   jupytext:
#     formats: py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: Python (ndmm-scrna)
#     language: python
#     name: ndmm-scrna
# ---

# %% [markdown]
# # 2025-07-02: Process BMMC Subjects [Stochastic]
# ### By [Aishwarya Chander](aishwarya.chander@alleninstitute.org), High Resolution Translational Immunology, Allen Institute for Immunology
#
# **Main aim**: \
# Plasma cells are highly hetrogeneous and donor specific. In this notebook, I create one object for each subject across all time points. I use this as the foundation to look into each subject for identifying clusters that may be tumor cells. 

# %% [markdown]
# ## 1. Import

# %%
import numpy as np
import pandas as pd
from scipy.stats import median_abs_deviation
import scanpy as sc
import os
import concurrent.futures

# %% [markdown]
# ## 2. Get data

# %%
adata = sc.read_h5ad('../../../data/rna/raw-files/all-bmmc-raw.h5ad')

# %%
bmmc_subject_dir = '../../../data/rna/bmmc-subjects/'
if not os.path.exists(bmmc_subject_dir):
    os.makedirs(bmmc_subject_dir)

# %% [markdown]
# ## 3. Process Subjects
# ### 3.1. Main Function

# %%
## The Median Absolute Deviation (MAD) is a robust measure of statistical dispersion.
## Compared to standard deviation, MAD better handles the zero-inflated, highly skewed nature of scRNA-seq data due to its robustness to outliers.
## It's calculated by:
# > Finding the median of a dataset
# > Computing absolute deviations from this median for each value
# > Taking the median of these absolute deviations

# > Stricter MAD thresholds (5) on full heterogeneous object; flags only, filtering deferred to downstream

from scipy.stats import median_abs_deviation
def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier


# %%
## Modified version of process_adata thats subject specific
def process_sub_adata(adata, subject_id, bmmc_subject_dir=bmmc_subject_dir):
  
    # Filter data for the specific subject
    sub_adata = adata[adata.obs['subject.subjectGuid'] == subject_id].copy()
    sub_adata.raw = sub_adata
    
    print(f"Processing subject {subject_id}")
    
    # Mitochondrial genes
    sub_adata.var["mt"] = sub_adata.var_names.str.startswith("MT-")
    # Ribosomal genes
    sub_adata.var["ribo"] = sub_adata.var_names.str.startswith(("RPS", "RPL"))
    # Hemoglobin genes
    sub_adata.var["hb"] = sub_adata.var_names.str.contains("^HB[^(P)]")
    
    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(sub_adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True)
    
    # Identify outliers
    sub_adata.obs["outlier"] = (
        is_outlier(sub_adata, "log1p_total_counts", 5)
        | is_outlier(sub_adata, "log1p_n_genes_by_counts", 5)
        | is_outlier(sub_adata, "pct_counts_in_top_20_genes", 5)
    )
    
    sub_adata.obs["mt_outlier"] = is_outlier(sub_adata, "pct_counts_mt", 3) | (
        sub_adata.obs["pct_counts_mt"] >= 6
    )
    
    # # Filter out doublets and outliers
    # sub_adata = sub_adata[sub_adata.obs["doublet_score"] < 0.3].copy()
    # sub_adata = sub_adata[(~sub_adata.obs.outlier) & (~sub_adata.obs.mt_outlier)].copy()
    
    # Store raw counts
    sub_adata.layers["counts"] = sub_adata.X.copy()
    
    # Normalizing to median total counts
    sc.pp.normalize_total(sub_adata)
    # Logarithmize the data
    sc.pp.log1p(sub_adata)
    
    # Identify highly variable genes
    sc.pp.highly_variable_genes(sub_adata, min_mean=0.0125, max_mean=3, min_disp=0.25)
    sub_adata = sub_adata[:, sub_adata.var_names[sub_adata.var['highly_variable']]].copy()
    
    # Scale the data
    sc.pp.scale(sub_adata, max_value=10, zero_center=False)
    # Perform PCA
    sc.tl.pca(sub_adata, svd_solver='arpack')
    sub_adata.obsm['X_pca_temp'] = sub_adata.obsm['X_pca']
    
    # Compute neighbors and UMAP
    sc.pp.neighbors(sub_adata, n_neighbors=50, use_rep='X_pca', n_pcs=20)
    sc.tl.umap(sub_adata, min_dist=0.45, random_state=0, n_components=2)
    sc.tl.leiden(sub_adata, resolution=1, n_iterations=2)
    
    # Save the processed data
    sub_adata.write_h5ad(f'{bmmc_subject_dir}bmmc-subject-{subject_id}.h5ad')
    print(f"Done saving subject {subject_id}")


# %% [markdown]
# ### 3.2. Run process loop
# #### Here, we set 'max_workers=5' to run 5 subjects in parallel as the files are now much larger than the base .h5s.

# %%
subject_ids = adata.obs['subject.subjectGuid'].unique()
with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
    futures = [executor.submit(process_sub_adata, adata, subject_id, bmmc_subject_dir=bmmc_subject_dir) for subject_id in subject_ids]
    concurrent.futures.wait(futures)
