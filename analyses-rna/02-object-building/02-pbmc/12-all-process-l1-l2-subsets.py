# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.2
#   kernelspec:
#     display_name: ndmm-scrna-envt
#     language: python
#     name: ndmm-scrna-envt
# ---

# %%
import pandas as pd
import numpy as np
import scanpy as sc
import scanpy.external as sce

sc.settings.n_jobs = 30
sc.settings.verbosity = 0

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)


# %%
import sys
sys.path.append("../../00-utilities/functions/python/")
from process_scrna_data import process_adata

# %% [markdown]
# ### Process L1 clusters with harmony

# %%
adata = sc.read_h5ad('../../../data/rna/final-objects/final-pbmc-raw.h5ad')

# %%
adata.obs['aifi_celltype_l1'].value_counts()

# %% jupyter={"outputs_hidden": true}
adata.raw = adata
clusters = list(adata.obs['aifi_celltype_l2'].value_counts().index)
clusters.reverse()

for cluster in clusters:
    print(f'Processing {cluster}')
    subset = adata[adata.obs['aifi_celltype_l2'] == cluster]
    subset = process_adata(subset, resolution=1, run_harmony=True,  run_rank_genes=True)
    subset.write(f'../../../data/rna/pbmc-celltypes/pbmc-{cluster}-processed-harmony.h5ad')
