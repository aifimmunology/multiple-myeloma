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
# This notebook partitions the non-plasma BMMC cells into major lineages (T cells, NK cells, B cells, monocytes, DCs, MSCs, progenitors) using marker gene expression across Leiden clusters, then independently reprocesses each lineage with Harmony batch correction. T and NK cells are further resolved by a second round of subsetting and reclustering to separate CD4 T, CD8 T, and NK populations, which are then saved as individual objects for downstream annotation.

# %% [markdown]
# ## 1. Imports

# %%
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

# %% [markdown]
# ## 2. Plot lineage markers

# %%
adata = sc.read_h5ad('../../../data/rna/bmmc/nonplasma-processed.h5ad')

# %% [markdown]
# ## 3. Annotate L1 cell types based on markers

# %%
sc.tl.dendrogram(adata, groupby='leiden')

dotplot = sc.pl.dotplot(
    adata,
    var_names=[
        # T cells
        "CD3D", "CD3E", "CD2", "TRAC", "IL7R", "CD4", "CD8A", "CCR5", "CXCR5", "DUSP1",

        # NK cells
        "SAMD1", "TBX21", "KLRB1", "NKG7", "GNLY", "KLRD1", "KLRF1", "FCGR3A", "PRF1",

        # B cells
        "CD19", "CD79A", "CD79B", "MS4A1", "CD22", "BANK1", "EBF1", "MZB1",

        # Monocytes
        "CD14", "LYZ", "FCGR3A", "S100A8", "S100A9", "VCAN", "FGR",

        # Dendritic cells (DCs)
        "CST3", "CLEC9A", "ITGAX", "BATF3", "CD1C", "FCER1A", "LILRA4", "AXL",

        # Mesenchymal stromal cells (MSCs)
        "LEPR", "ENG", "THY1", "PDGFRA", "NT5E", "CXCL12",

        # Progenitors
        "PTPRC", "MKI67", "CD34", "GATA2", "KIT", "SOX4", "TOP2A", "HES1",
    ],
    groupby='leiden',
    standard_scale='var',
    cmap='Reds',
    dendrogram=True,
    show=False 
)

# Save the figure
plt.savefig("figures/dotplot_celltype_markers.png", dpi=300, bbox_inches="tight")
plt.close()

# %%
leiden_to_label = {'0': 't_cells',
 '1': 't_cells',
 '2': 't_cells',
 '3': 't_cells',
 '4': 't_cells',
 '5': 't_cells',
 '6': 't_cells',
 '7': 't_cells',
 '8': 'monocytes',
 '9': 'b_cells',
 '10': 'b_cells',
 '11': 'b_cells',
 '12': 't_cells',
 '13': 't_cells',
 '14': 't_cells',
 '15': 'nk_cells',
 '16': 'nk_cells',
 '17': 'monocytes',
 '18': 'progenitors',
 '19': 'dcs',
 '20': 'monocytes',
 '21': 't_cells',
 '22': 'mscs',
 '23': 'dcs',
 '24': 'monocytes',
 '25': 't_cells',
 '26': 't_cells',
 '27': 't_cells',
 '28': 't_cells',
 '29': 't_cells',
 '30': 't_cells',
 '31': 't_cells',
 '32': 't_cells',
 '33': 't_cells',
 '34': 'nk_cells',
 '35': 'nk_cells',
 '36': 'nk_cells',
 '37': 'b_cells',
 '38': 'b_cells',
 '39': 'progenitors',
 '40': 'progenitors',
 '41': 'progenitors',
 '42': 'b_cells',
 '43': 'b_cells',
 '44': 'progenitors',
 '45': 't_cells',
 '46': 'nk_cells',
 '47': 'b_cells'}

adata.obs['temp.aifi_l1'] = adata.obs['leiden'].map(leiden_to_label).astype('category').cat.remove_unused_categories()

# %%
sc.pl.umap(
    adata,
    color=['temp.aifi_l1'],
    size=0.5,
    legend_loc='on data',
    show=False,
    save="_partition_l1.png"
)

# %%
celltypes = list(adata.obs['temp.aifi_l1'].value_counts().index)
celltypes.reverse()

for celltype in celltypes:
    print(f'Processing {celltype}')
    subset = adata[adata.obs['temp.aifi_l1'] == celltype]

    if celltype == 'b_cells':
        subset = subset.raw.to_adata()   
        igl_genes = [gene for gene in subset.var_names if gene.startswith('IGL')]
        igk_genes = [gene for gene in subset.var_names if gene.startswith('IGK')]
        ighc_genes = [gene for gene in subset.var_names if gene.startswith('IGH')]
        all_ig_genes = igl_genes + igk_genes + ighc_genes
        filtered_genes = [gene for gene in subset.var_names if gene not in all_ig_genes]
        subset = subset[:, filtered_genes]
        subset.raw = subset

    subset = process_adata(subset, resolution=1, run_harmony=True)
    subset.write(f'../../../data/rna/bmmc/{celltype}-processed-v1.h5ad')

# %% [markdown]
# ### 3.1. Subset T and NK cells and pull apart locally
# > This needs to be done as the lines between CD8+ T cells and NK Cells are hard to draw in the total object. The added local resolution exascerbates the differences between T and NK cells.

# %%
sub_t_nk = adata[adata.obs['temp.aifi_l1'].isin(['nk_cells', 't_cells'])].copy()
sub_t_nk = process_adata(sub_t_nk, resolution=1.5, run_harmony=True)
sub_t_nk.write(f'../../../data/rna/bmmc/t_nk_cells-processed-v1.h5ad')

# %%
sub_t_nk = sc.read_h5ad('../../../data/rna/bmmc/t_nk_cells-processed-v1.h5ad')

# %%
leiden_to_label ={'0': 't_cd8',
 '1': 't_cd8',
 '2': 't_cd8',
 '3': 'nk_cells',
 '4': 'nk_cells',
 '5': 't_cd4',
 '6': 't_cd4',
 '7': 't_cd4',
 '8': 't_cd4',
 '9': 't_cd4',
 '10': 't_cd4',
 '11': 't_cd8',
 '12': 't_cd4',
 '13': 't_cd4',
 '14': 't_cd4',
 '15': 't_cd4',
 '16': 't_cd4',
 '17': 't_cd8',
 '18': 't_cd4',
 '19': 't_cd4',
 '20': 'nk_cells',
 '21': 'nk_cells'}

sub_t_nk.obs['partition'] = sub_t_nk.obs['leiden'].map(leiden_to_label).astype('category').cat.remove_unused_categories()

# %%
celltypes = list(sub_t_nk.obs['partition'].value_counts().index)
celltypes.reverse()

for celltype in celltypes:
    print(f'Processing {celltype}')
    subset = sub_t_nk[sub_t_nk.obs['partition'] == celltype]
    subset = process_adata(subset, resolution=1.5, run_harmony=True)
    subset.write(f'../../../data/rna/bmmc/{celltype}-processed-v1.h5ad')
