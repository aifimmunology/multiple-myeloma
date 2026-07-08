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
# # 2025-07-01: Process BMMC object [Stochastic]
# ### By [Aishwarya Chander](aishwarya.chander@alleninstitute.org), High Resolution Translational Immunology, Allen Institute for Immunology
#
# **Main aim**: \
# Here, I'll process my object prior to any QC or filtering to get a glimpse of what my data looks like. This allows me to craft a plan thats specific to the issues in this dataset for further steps. The `process_adata` function is common across all future data processing. Setting `sc.settings.n_jobs = 30` is a way to increase or limit the default number of jobs or CPUs for parallel computing in all Scanpy Processes. 

# %% [markdown]
# ## 1. Imports

# %% editable=true slideshow={"slide_type": ""}
import os
import scanpy as sc
import scanpy.external as sce
import matplotlib.pyplot as plt
import numpy as np

sc.settings.n_jobs = 30
sc.settings.verbosity = 0

plt.rcParams['figure.dpi'] = 80
plt.rcParams['figure.figsize'] = (5, 5)

import sys
sys.path.insert(0, '../../00-utilities/functions/python/')
from process_scrna_data import process_adata

# %%
plasma_dir = '../../../data/rna/plasma/'
if not os.path.exists(plasma_dir):
    os.makedirs(plasma_dir)

bmmc_dir = '../../../data/rna/bmmc/'
if not os.path.exists(bmmc_dir):
    os.makedirs(bmmc_dir)

# %% [markdown]
# ## 2. Process all BMMCs

# %%
adata = sc.read_h5ad('../../../data/rna/raw-files/all-bmmc-raw.h5ad')
adata.raw = adata

# %%
adata = process_adata(adata)

# %% [markdown]
# ## 3. Make some plots

# %% [markdown]
# ### 3.1. View plasma markers 

# %%
sc.pl.umap(
    adata,
    color=['SDC1', 'JUN', 'TNFRSF17', 'TXNDC5', 'SLAMF7', 'PRDM1', 'MZB1'],
    cmap='Reds',
    size=0.5,
    vmax='p99',
    use_raw=False,
    show=False,
    save="_plasma_genes.png" 
)

# %% [markdown]
# ### 3.2. View leiden clusters

# %%
sc.pl.umap(
    adata,
    color=['leiden'],
    size=0.5,
    legend_loc='on data',
    show=False,
    save="_leiden_clusters.png" 
)

# %% [markdown]
# ## 4. Annotate potential plasma cells
# #### (Stochastic)

# %%
#### Based on the saved plasama gene expression and leiden cluster UMAPs, we can now annotate a group of cells as tumor like. This needs to be adjuster further as it is currently imprecise. 
adata.obs['possibly_tumor'] = np.where(
    adata.obs['leiden'].isin(['19', '26', '34', '35']),
    'maybe plasma',
    'non-plasma'
)

# %% [markdown]
# ## 5. View relevant features

# %%
sc.pl.umap(
    adata,
    color=['healthy_l1', 'manual.category', 'possibly_tumor'],
    size=0.5,
    legend_loc='on data',
    show=False,
    save="_categories.png"
)

# %% [markdown]
# ## 6. Save object and metadata

# %%
adata.write(bmmc_dir+'all-bmmc-processed.h5ad')
adata.obs.to_parquet(bmmc_dir+'all-bmmc-processed-metadata.parquet')

# %% [markdown]
# ## 7. Save preliminary plasma cell UUIDs

# %%
subset_tumor = adata.obs.loc[adata.obs['possibly_tumor'] == 'maybe plasma', 'cell_uuid']
subset_tumor.to_frame(name="possibly_tumor").to_parquet(plasma_dir+'preliminary_plasma_uuids.parquet')
