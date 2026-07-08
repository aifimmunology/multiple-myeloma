# ---
# jupyter:
#   jupytext:
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
# %%
# adata = sc.read_h5ad('../../../data/rna/final-objects/final-pbmc-raw.h5ad')
mask = adata.obs['aifi_celltype_l1'].isin(['mono', 'dc'])
adata = adata[mask, :]
adata = process_adata(adata, resolution=1, run_harmony=True,  run_rank_genes=True)
adata.write(f'../../../data/rna/pbmc-subsets/pbmc-mono_dc-processed-harmony.h5ad')

# %%
## Format for scarf
# adata = sc.read_h5ad('../../../data/rna/final-objects/final-pbmc-raw.h5ad')
# adata.var["gene_id"] = adata.var.index.astype(str)
# adata.write_h5ad("../../../data/rna/final-objects/final-pbmc-raw-scarf.h5ad")

## Extract only T cells
# adata = sc.read_h5ad('../../../data/rna/final-objects/final-pbmc-raw.h5ad')
# adata = adata[adata.obs['aifi_celltype_l1'] == 't_cell']
# adata.write('../../../data/rna/pbmc-subsets/all-pbmc-t-cells.h5ad')
# sc.pp.subsample(adata, n_obs=50000, random_state=0)
# adata.write('../../../data/rna/pbmc-subsets/all-pbmc-t-cells-50k.h5ad')