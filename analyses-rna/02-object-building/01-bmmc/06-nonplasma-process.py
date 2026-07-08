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
# # 2025-07-04: Build Non-plasma cells Object [Stochastic]
# ### By [Aishwarya Chander](aishwarya.chander@alleninstitute.org), High Resolution Translational Immunology, Allen Institute for Immunology
# **Main aim**: 
# After removing plasma cells and predicted doublets from the full BMMC object, this notebook performs QC filtering on the remaining non-plasma cells using MAD-based outlier detection (total counts, gene counts, top gene fraction) and a mitochondrial percentage threshold. The filtered object is then processed through the standard scanpy pipeline (normalization, HVG selection, PCA, UMAP, Leiden clustering at resolution 2) to generate a clean, clustered non-plasma BMMC object for downstream annotation.

# %% [markdown]
# ## 1. Imports

# %% editable=true slideshow={"slide_type": ""}
import anndata
import scanpy as sc
import matplotlib.pyplot as plt
import os
import scanpy.external as sce
import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import median_abs_deviation

sc.settings.n_jobs = 30
sc.settings.verbosity = 0

plt.rcParams['figure.dpi'] = 80
plt.rcParams['figure.figsize'] = (5, 5)

# %%
import sys
sys.path.insert(0, '../../00-utilities/functions/python/')
from process_scrna_data import process_adata

# %%
# The Median Absolute Deviation (MAD) is a robust measure of statistical dispersion.
# Compared to standard deviation, MAD better handles the zero-inflated, highly skewed nature of scRNA-seq data due to its robustness to outliers.
# It's calculated by:
# > Finding the median of a dataset
# > Computing absolute deviations from this median for each value
# > Taking the median of these absolute deviations

# > # Looser MAD thresholds (7) appropriate for non-plasma subset where distributions are tighter post-cleanup

def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier


# %% [markdown]
# ## 2. Process object
# > Start with the raw data object. Then, add in the plasma cell IDs from the plasma cell analysis. 

# %%
adata = sc.read_h5ad('../../../data/rna/raw-files/all-bmmc-raw.h5ad')
adata.raw = adata

# %%
plasma_all = pd.read_parquet('../../../data/rna/plasma/all-plasma-metadata.parquet')
subset_tumor = pd.read_parquet('../../../data/rna/plasma/preliminary_plasma_uuids.parquet')


drop_plasma_uuids = set(plasma_all['cell_uuid']) | set(subset_tumor['possibly_tumor'])

# %%
adata = adata[~adata.obs['cell_uuid'].isin(list(drop_plasma_uuids))].copy()

# %% [markdown]
# ## 3. QC and filtering

# %%
# mitochondrial genes
adata.var["mt"] = adata.var_names.str.startswith("MT-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes.
adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]"))

sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
)

# %%
adata.obs["outlier"] = (
    is_outlier(adata, "log1p_total_counts", 7)
    | is_outlier(adata, "log1p_n_genes_by_counts", 7)
    | is_outlier(adata, "pct_counts_in_top_20_genes", 7)
)
adata.obs.outlier.value_counts()

# %%
adata.obs["mt_outlier"] = adata.obs["pct_counts_mt"] >  6
adata.obs.mt_outlier.value_counts()

# %%
sum((adata.obs.outlier) & (adata.obs.mt_outlier))

# %%
print(f"Total number of cells: {adata.n_obs}")

adata = adata[adata.obs["doublet_score"] <= 0.35].copy()
adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()

print(f"Number of cells after filtering of low quality cells: {adata.n_obs}")

# %% [markdown]
# ## 4. Processing

# %%
adata = process_adata(adata, resolution=2)

# %% [markdown]
# ### Save object and metadata

# %%
adata.write('../../../data/rna/bmmc/nonplasma-processed.h5ad')
adata.obs.to_parquet('../../../data/rna/bmmc/nonplasma-metadata.parquet')
