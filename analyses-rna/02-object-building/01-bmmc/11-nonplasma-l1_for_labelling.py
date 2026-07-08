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
#     display_name: ndmm-scrna-envt
#     language: python
#     name: ndmm-scrna-envt
# ---

# %% [markdown]
# # 2025-07-27: Reprocess BMMC cell types at L1 and L2 for manual label review [Stochastic]
# ### By [Aishwarya Chander](aishwarya.chander@alleninstitute.org), High Resolution Translational Immunology, Allen Institute for Immunology
# **Main aim**: This notebook takes the refined BMMC object (with KNN-imputed labels from step 10) and reprocesses each L2 and L1 cell type subset independently with Harmony batch correction and Leiden clustering, generating per-lineage objects for manual inspection and label cleanup. Cell types with fewer than 100 cells are skipped.

# %% [markdown]
# ## 1. Imports

# %%
import glob
import os

import anndata
import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import concurrent.futures

sc.settings.n_jobs = 30
sc.settings.verbosity = 0

plt.rcParams["figure.dpi"] = 100
plt.rcParams["figure.figsize"] = (6, 6)

# %%
import sys
sys.path.insert(0, '../../00-utilities/functions/python/')
from process_scrna_data import process_adata

# %% [markdown]
# ## 2. Process object

# %%
adata = sc.read_h5ad('../../../data/rna/bmmc-outliers/bmmc-refined-labels.h5ad')

# %% [markdown]
# ### 2.1. Process at L2

# %%
clusters = list(adata.obs['aifi_celltype_l2_knn'].value_counts(ascending=True).index)

for cluster in clusters:
    subset = adata[adata.obs['aifi_celltype_l2_knn'] == cluster]
    
    if subset.n_obs < 100:
        print(f'Skipping {cluster} (n={subset.n_obs}) — too few cells')
        continue

    print(f'Processing {cluster}')
    subset = process_adata(subset, resolution=1, run_harmony=True, run_rank_genes=True)
    subset.write(f'../../../data/rna/bmmc-cleanup/cleanup-l2-{cluster}.h5ad')

# %% [markdown]
# ### 2.2. Process at L2

# %%
clusters = list(adata.obs['aifi_celltype_l1_knn'].value_counts(ascending=True).index)

for cluster in clusters:
    subset = adata[adata.obs['aifi_celltype_l1_knn'] == cluster]
    
    if subset.n_obs < 100:
        print(f'Skipping {cluster} (n={subset.n_obs}) — too few cells')
        continue

    print(f'Processing {cluster}')
    subset = process_adata(subset, resolution=1, run_harmony=True, run_rank_genes=True)
    subset.write(f'../../../data/rna/bmmc-cleanup/cleanup-l1-{cluster}.h5ad')
