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

# %%
import pandas as pd
import hisepy as hp
import numpy as np
import scanpy as sc
import scanpy.external as sce
import matplotlib.pyplot as plt
sc.settings.n_jobs = 60


# %%
import sys
sys.path.append("../../00-utilities/functions/python/")
from process_scrna_data import process_adata

# %%
adata = sc.read_h5ad('../../../data/rna/pbmc-subsets/all-pbmc-filtered-500k.h5ad')
adata.raw = adata

# %%
adata = process_adata(
    adata,
    resolution=0.1,
)

# %%
adata.write('../../../data/rna/pbmc/adata-pbmc-500k-umap.h5ad')

# %%
adata

# %%
l1_cmap = {
    'T cell': '#5480A3',
    'B cell': '#A3787B',
    'NK cell': '#FF7F55',
    'Monocyte': '#B75228',
    'DC': '#F8C755',
    'ILC': '#F4AD3A',
    'Progenitor cell': '#A31001',
    'Platelet': '#5B3930',
    'Erythrocyte': '#1F1515',
    'Other': '#1F1515'
}

l2_cmap = {
    'Naive CD4 T cell': '#597FC6',
    'Memory CD4 T cell': '#233C5B',
    'Treg': '#001A9A',
    'DN T cell': '#656E68',
    'Proliferating T cell': '#AEAEB0',
    'Naive CD8 T cell': '#63A686',
    'Memory CD8 T cell': '#314937',
    'CD8aa': '#7ECDAD',
    'gdT': '#72AF5A',
    'MAIT': '#1E8F64',
    'Transitional B cell': '#9C6469',
    'Naive B cell': '#D5ABAB',
    'Memory B cell': '#B16B72',
    'Effector B cell': '#F4CEC3',
    'Plasma cell': '#784D47',
    'CD56bright NK cell': '#FF7F55',
    'CD56dim NK cell': '#E16040',
    'Proliferating NK cell': '#CC3B18',
    'CD14 monocyte': '#F6831B',
    'Intermediate monocyte': '#B75228',
    'CD16 monocyte': '#FCAF87',
    'cDC2': '#BE6E23',
    'cDC1': '#F8C755',
    'pDC': '#FBB64E',
    'ASDC': '#8A4E1C',
    'Other_DC': '#DDAA66',
    'ILC': '#F4AD3A',
    'Progenitor cell': '#A31001',
    'Platelet': '#5B3930',
    'Erythrocyte': '#1F1515'
}

l3_cmap = {
    'SOX4+ naive CD4 T cell': '#2B5777',
    'Core naive CD4 T cell': '#597FC6',
    'ISG+ naive CD4 T cell': '#072E77',
    'CM CD4 T cell': '#6793A4',
    'GZMB- CD27+ EM CD4 T cell': '#406B9C',
    'GZMB- CD27- EM CD4 T cell': '#233C5B',
    'ISG+ memory CD4 T cell': '#5480A3',
    'KLRF1- GZMB+ CD27- memory CD4 T cell': '#628088',
    'Naive CD4 Treg': '#8EB5D7',
    'Memory CD4 Treg': '#001A9A',
    'KLRB1+ memory CD4 Treg': '#3B76C1',
    'GZMK+ memory CD4 Treg': '#A3B7C2',
    'Memory CD8 Treg': '#C3D1DA',
    'KLRB1+ memory CD8 Treg': '#C7E0C6',
    'DN T cell': '#656E68',
    'Proliferating T cell': '#AEAEB0',
    'SOX4+ naive CD8 T cell': '#465F32',
    'Core naive CD8 T cell': '#63A686',
    'ISG+ naive CD8 T cell': '#869A94',
    'CM CD8 T cell': '#6D9F5E',
    'GZMK+ CD27+ EM CD8 T cell': '#95BDA1',
    'KLRF1- GZMB+ CD27- EM CD8 T cell': '#59895F',
    'GZMK- CD27+ EM CD8 T cell': '#314937',
    'KLRF1+ GZMB+ CD27- EM CD8 T cell': '#4E695B',
    'ISG+ memory CD8 T cell': '#94C0AD',
    'CD8aa': '#7ECDAD',
    'SOX4+ Vd1 gdT': '#B6D7C8',
    'Naive Vd1 gdT': '#537E41',
    'KLRF1- effector Vd1 gdT': '#A39E78',
    'KLRF1+ effector Vd1 gdT': '#7CB38E',
    'GZMK+ Vd2 gdT': '#72AF5A',
    'GZMB+ Vd2 gdT': '#AECFB1',
    'CD8 MAIT': '#1E8F64',
    'CD4 MAIT': '#6F723F',
    'ISG+ MAIT': '#BBB788',
    'Transitional B cell': '#9C6469',
    'Core naive B cell': '#D5ABAB',
    'ISG+ naive B cell': '#D9A996',
    'Early memory B cell': '#362027',
    'Core memory B cell': '#B16B72',
    'Activated memory B cell': '#60545A',
    'Type 2 polarized memory B cell': '#4E2927',
    'CD95 memory B cell': '#EBC3A7',
    'CD27+ effector B cell': '#F4CEC3',
    'CD27- effector B cell': '#A3787B',
    'Plasma cell': '#784D47',
    'CD56bright NK cell': '#FF7F55',
    'GZMK+ CD56dim NK cell': '#A31001',
    'GZMK- CD56dim NK cell': '#E16040',
    'Adaptive NK cell': '#781B15',
    'ISG+ CD56dim NK cell': '#EEA172',
    'Proliferating NK cell': '#CC3B18',
    'Core CD14 monocyte': '#F6831B',
    'ISG+ CD14 monocyte': '#DE5A07',
    'IL1B+ CD14 monocyte': '#BE4C26',
    'Intermediate monocyte': '#B75228',
    'Core CD16 monocyte': '#FCAF87',
    'ISG+ CD16 monocyte': '#DB6A2E',
    'C1Q+ CD16 monocyte': '#F4BA81',
    'CD14+ cDC2': '#FEDC7A',
    'HLA-DRhi cDC2': '#BE6E23',
    'ISG+ cDC2': '#8C6527',
    'cDC1': '#F8C755',
    'pDC': '#FBB64E',
    'ASDC': '#8A4E1C',
    'Other_DC': '#DDAA66',
    'ILC': '#F4AD3A',
    'CLP cell': '#FCDB97',
    'CMP cell': '#5B3006',
    'BaEoMaP cell': '#3D1F15',
    'Platelet': '#5B3930',
    'Erythrocyte': '#1F1515'
}

# %%
plt.rcParams['figure.dpi'] = 120
plt.rcParams['figure.figsize'] = (8, 8)

# %%
sc.pl.umap(
    adata,
    color=['aifi_label_l2'], 
    frameon='False',
    palette=l2_cmap,
    size=1
)

# %%
sc.pl.umap(
    adata,
    color=['aifi_label_l2'], 
    frameon='False',
    palette=l2_cmap,
    size=1
)

# %%
# CMV serostatus per subject (source: clinical metadata)
cmv_df = pd.read_csv('../../../data/rna/metadata/subject_cmv_status.csv')
cmv_dict = dict(zip(cmv_df['subject_id'], cmv_df['cmv_status']))
adata.obs["subject.cmv"] = adata.obs["subject.subjectGuid"].map(cmv_dict)

# %%
sc.pl.umap(
    adata,
    color=['subject.cmv'],
    frameon=False,
    size=1, save='_pbmc_cmv.png' )

# %%
sc.pl.umap(
    adata,
    color=['leiden'],
    size=1,
    legend_loc='on data',
    show=False,
    save='_pbmc_leiden_clusters.png' 
)

# %%
sc.pl.umap(
    adata,
    color=['aifi_label_l1'],
    palette=l1_cmap,
    size=1,
    legend_loc='on data',
    show=False,
    save='_pbmc_leiden_l1.png' 
)

# %%
sc.pl.umap(
    adata,
    color=['aifi_label_l1'],
    palette=l1_cmap,
    size=1,
    show=False,
    save='_pbmc_leiden_l1-legend.png' 
)

# %%
sc.pl.umap(
    adata,
    color=['aifi_label_l2'],
    palette=l2_cmap,
    size=1,
    legend_loc='on data',
    show=False,
    save='_pbmc_leiden_l2.png' 
)

# %%
sc.pl.umap(
    adata,
    color=['aifi_label_l2'],
    palette=l2_cmap,
    size=1,
    show=False,
    save='_pbmc_leiden_l2-legend.png' 
)

# %%
sc.pl.umap(
    adata,
    color=['aifi_label_l3'],
    palette=l3_cmap,
    size=1,
    legend_loc='on data',
    show=False,
    save='_pbmc_leiden_l3.png' 
)

# %%
sc.pl.umap(
    adata,
    color=['aifi_label_l3'],
    palette=l3_cmap,
    size=1,
    show=False,
    save='_pbmc_leiden_l3-legend.png' 
)

# %%
sc.pl.umap(
    adata,
    color=['cohort.cohortGuid'],
    size=1,
    show=False,
    save='_pbmc_cohort-legend.png' 
)
