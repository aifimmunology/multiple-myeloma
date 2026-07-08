# ---
# jupyter:
#   jupytext:
#     formats: py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.3
#   kernelspec:
#     display_name: Python (scNMF)
#     language: python
#     name: scnmf
# ---

# %% [markdown]
# # 2025-09-08: PBMC - Run scNMF on CD4 T cells - MM Cells
# ### By [Aishwarya Chander](aishwarya.chander@alleninstitute.org), High Resolution Translational Immunology, Allen Institute for Immunology
#
# **Main aim**: \
# In this notebook, I will use the method described by [Yasumizu et al.](https://github.com/yyoshiaki/NMFprojection) to perform NMF decomposition with a pre-computed matrix W from Human CD4 T cells to identify niche T cell states within all CD4 + T cells in the PBMC dataset. This notebook uses the kernel defined in the `rna-setup/ac-nfm-envt.yml` file. 
# > Original tutorial by [Pravina Venkatesan](https://github.com/aifimmunology/ALTRA2_Analysis/blob/primary_analysis_PV/scRNA/NMF_Analysis/00_NMF_Projection.ipynb)
#
# NMF is a memory intensive step. It is recommended to run this script outside an `.ipynb` notebook. 

# %% [markdown]
# ### 1. Imports

# %%
import scanpy as sc
import os
import pandas as pd
import numpy as np
from scipy.io import mmread
import scipy.sparse as sp
import scanpy as sc
from NMFproj import *
import gc

import sys
sys.path.append("../../00-utilities/functions/python/")
from process_scrna_data import process_adata

# %% [markdown]
# ### 2. Subset all T cells to CD4 T cells 

# %%
adata = sc.read_h5ad('../../../data/rna/pbmc-subsets/all-pbmc-t-cells.h5ad')

t_cd4_pos = ['t_reg_cd4_naive',
             't_reg_cd4_memory-klrb1.pos',
             't_reg_cd4_memory',
             't_cd4_memory_central',
             't_cd4_memory_effector-cd27.neg-gzmb.neg',
             't_cd4_naive-core',
             't_cd4_memory_effector-cd27.pos-gzmb.neg',
             't_cd4_naive-isg.pos',
             't_cd4_memory-isg.pos',
             't_mait-cd4',
             't_cd4_naive-sox4.pos']

adata = adata[adata.obs['aifi_celltype_l3'].isin(t_cd4_pos)]
adata = adata[adata.obs['manual.category']
              == 'tumor_pbmc']  # Keep only MM T cells

adata.raw = adata  # add in raw adata

adata.write('../../../data/rna/pbmc-subsets/pbmc-t-cd4-mm-input-nmf.h5ad')

# Scanpy Normalize to avoid sum step using the matrix
sc.pp.normalize_total(adata, target_sum=1e4, copy=False)
sc.pp.log1p(adata, copy=False)

# %%
# Grab normalized matrix as input for NMF
X = adata.to_df().T
print('Raw RNA inputs defined. Object Shape:', X.shape)

# %%
# Save some memory before moving forward
adata = None
gc.collect()

# %% [markdown]
# ### 3. Run NMF projection steps using precomputed weights

# %%
## Get the Pre-computed CD4 T cell weight matrix from paper (github - https://github.com/yyoshiaki/NMFprojection/blob/main/data/NMF.W.CD4T.csv.gz)
fixed_W = pd.read_csv("NMF_W_CD4T.csv")
fixed_W = fixed_W.rename(columns={'Unnamed: 0': ''})
fixed_W.set_index(fixed_W.columns[0], inplace=True)

# %%
## Run the NMFproj function; here, I'm using the raw counts so they're not normalized
X_norm, X_trunc, df_H, fixed_W_trunc = NMFproj(X, fixed_W, return_truncated=True, normalized=True)
print('NMF projection computed.')

# %% [markdown]
# ### 4. Add in computed projections into adata

# %%
adata = sc.read_h5ad('../../../data/rna/pbmc-subsets/pbmc-t-cd4-mm-input-nmf.h5ad')

# %%
df_ev = calc_EV(X_trunc, fixed_W_trunc, df_H)
df_H.iloc[:, :2]

# %%
index_mapping = {
    'NMF_0': 'NMF0_Cytotoxic',
    'NMF_1': 'NMF1_Treg',
    'NMF_2': 'NMF2_Th17',
    'NMF_3': 'NMF3_Naive',
    'NMF_4': 'NMF4_Act',
    'NMF_5': 'NMF5_Th2',
    'NMF_6': 'NMF6_Tfh',
    'NMF_7': 'NMF7_IFN',
    'NMF_8': 'NMF8_Cent_Mem',
    'NMF_9': 'NMF9_Thymic_Emi',
    'NMF_10': 'NMF10_Tissue',
    'NMF_11': 'NMF11_Th1'
}
## Normalize NMF feature matrix and add to adata.obs
df_H.index = df_H.index.map(index_mapping)
df_H_norm = (df_H.T / df_H.max(axis=1))
adata.obs = pd.merge(adata.obs, df_H_norm, left_index=True, right_index=True)

# %%
## List of the NMF components
nmf_components = list(index_mapping.values())

## Add in the NMF scores as a matrix for each cell
combined_matrix = adata.obs[nmf_components].values
combined_matrix
adata.obsm['nmf_factors'] = combined_matrix

## NMF is projected with a range of values for each cell, I'll label each cell with the NMF component that has the highest score
nmf_values = adata.obs[nmf_components].values
max_component_indices = np.argmax(nmf_values, axis=1)
adata.obs['nmf_cell_labels'] = [nmf_components[i] for i in max_component_indices]

# %%
adata.write('../../../data/rna/pbmc-subsets/pbmc-t-cd4-mm-nmf-raw.h5ad')
adata.obs.to_parquet('../../../data/rna/pbmc-subsets/pbmc-t-cd4-mm-nmf-metadata.parquet')
print('Raw files saved.')

# %% [markdown]
# ### 5. Process object on the PCA space

# %%
adata = process_adata(adata, run_harmony=True, resolution=0.5)
adata.obsm['base_umap'] = adata.obsm['X_umap'].copy()
adata.obsm['base_tsne'] = adata.obsm['X_tsne'].copy()
print('Regular processing done.')

# %% [markdown]
# ### 6. Process object on the NMF space

# %%
sc.pp.neighbors(adata, n_neighbors=15, use_rep='nmf_factors', key_added='nmf_neighbors')
sc.tl.umap(adata, neighbors_key='nmf_neighbors')
adata.obsm['nmf_umap'] = adata.obsm['X_umap'].copy()

sc.tl.tsne(adata, use_rep='nmf_factors')
adata.obsm['nmf_tsne'] = adata.obsm['X_tsne'].copy()

sc.tl.leiden(adata, resolution=1, key_added='nmf_leiden',
             n_iterations=2, neighbors_key='nmf_neighbors')

print('NMF processing done.')

# %%
## Save processed object
adata.write('../../../data/rna/pbmc-subsets/pbmc-t-cd4-mm-nmf-processed.h5ad')
print('Processed files saved.')
