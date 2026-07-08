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
# # 2025-07-10: CellTypist Annotations
# ### By [Aishwarya Chander](aishwarya.chander@alleninstitute.org), High Resolution Translational Immunology, Allen Institute for Immunology
#
# **Main aim**: This includes 3 sets of CellTypist models; an external healthy BMMC dataset, internal healthy BMMC dataset and AIFI reference PBMC dataset.

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## 1. Imports and functions

# %% editable=true slideshow={"slide_type": ""}
import hisepy as hp

import os
import gc
import glob
import shutil

import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as scs
import anndata

from multiprocessing import Pool

import celltypist
from celltypist import models

models.download_models(
    model=['Immune_All_High.pkl',
           'Immune_All_Low.pkl']
)

## Object builders
import sys
sys.path.insert(0, '../../00-utilities/functions/python/')

from h5_data_reader import read_h5_anndata

# %%
output_base_path = '../../../data/rna/celltypist/labels/'
if not os.path.exists(output_base_path):
    os.makedirs(output_base_path)

model_base_path = '../../../data/rna/celltypist/models/'
if not os.path.exists(model_base_path):
    os.makedirs(model_base_path)

# %% [markdown]
# ## 2. Get metadata
#
# Pull metadata from HISE to grab all relevant file paths for label prediction. These already have relevant H5AD file paths in them relative to this notebooks path.

# %% editable=true slideshow={"slide_type": ""}
metadata = pd.read_csv('../../../data/rna/metadata/scrna_metadata.csv')

# %% [markdown]
# ## 3. BMMC Celltypst Labelling
#
# First, we'll download our 3 reference CellTypist models into our directory. Then, we'll parallelize the celltype annotation across all our samples.

# %%
# Model setup
## Download CellTypist reference models from HISE and place in models directory
import shutil

model_base_path = '../../../data/rna/celltypist/models/'
os.makedirs(model_base_path, exist_ok=True)

# BMMC references (Triana et al. external + AIFI internal healthy)
bmmc_models = {
    'bmmc_ext_reference_l1.pkl': '6b7955ce-d0d5-4ad8-8d2a-bfe0c8f694ad',
    'bmmc_ext_reference_l2.pkl': '634b9fd3-596b-47c4-9d2e-77e238d54ffd',
    'bmmc_ext_reference_l3.pkl': 'b9dddba8-6fc8-4986-a4ac-c64efbe40a3f',
    'bmmc_healthy_reference_l1.pkl': 'bd7aa60a-1df3-4252-a349-482c80ec4d7a',
    'bmmc_healthy_reference_l2.pkl': '384a7518-a075-4707-827e-a8251c17f8b0',
    'bmmc_healthy_reference_l3.pkl': 'ff0a57c9-3dc5-455f-9d54-a5c43bad26c9',
}

# AIFI PBMC references
pbmc_models = {
    'ref_pbmc_clean_celltypist_model_flex-features_AIFI_L1_2024-04-18.pkl': '05ee02f5-baa3-44dc-a712-40d4075c1cfb',
    'ref_pbmc_clean_celltypist_model_flex-features_AIFI_L2_2024-04-19.pkl': '3fe50df2-95a3-4bd8-8822-b338f19e412c',
    'ref_pbmc_clean_celltypist_model_flex-features_AIFI_L3_2024-04-19.pkl': 'ee1e4556-0768-4eb6-9f49-c6b9567df331',
}

all_models = {**bmmc_models, **pbmc_models}

# Download from HISE and move to model directory
cached_paths = hp.cache_files(list(all_models.values()))

for model_name, uuid in all_models.items():
    dest = os.path.join(model_base_path, model_name)
    if not os.path.exists(dest):
        # Find the downloaded file matching this UUID
        src = [p for p in cached_paths if uuid in p]
        if src:
            shutil.copy2(src[0], dest)

print(f"Models available in {model_base_path}: {os.listdir(model_base_path)}")


# %% [markdown]
# ### Loops through all 3 loops and collects celltypist predictions

# %% editable=true slideshow={"slide_type": ""}
def celltypist_processing_bmmc(args):
    # Read data
    data_file, model_base_path, output_base_path = args
    pbmc = read_h5_anndata(data_file)

    # Processing testing data
    sample_id = pbmc.obs['pbmc_sample_id'].unique().tolist()[0]
    pbmc.obs.index = pbmc.obs['barcodes']

    # Normalization and log transformation
    sc.pp.normalize_total(pbmc, target_sum=1e4)
    sc.pp.log1p(pbmc)

    # Annotations
    codes = ['ext', 'healthy']
    levels = ['l1', 'l2', 'l3']
    predictions = {}

    for code in codes:
        for level in levels:
            model_file = f'{model_base_path}bmmc_{code}_reference_{level}.pkl'
            predictions[f'{level}_{code}'] = celltypist.annotate(
                pbmc, model=model_file)

            # Write out labels
            predictions[f'{level}_{code}'].predicted_labels.reset_index().to_csv(
                f'{output_base_path}{sample_id}_{level}_{code}_predicted_labels.csv')
            # Write out probability matrix
            predictions[f'{level}_{code}'].probability_matrix.reset_index().to_parquet(
                f'{output_base_path}{sample_id}_{level}_{code}_probability_matrix.parquet')
            # Write out decision matrix
            predictions[f'{level}_{code}'].decision_matrix.reset_index().to_parquet(
                f'{output_base_path}{sample_id}_{level}_{code}_decision_matrix.parquet')


# %% jupyter={"outputs_hidden": true}
args_list = [(file, model_base_path, output_base_path)
             for file in metadata['manual.file_paths']]

with Pool(processes=30) as pool:
    pool.map(celltypist_processing_bmmc, args_list)


# %% [markdown]
# ## 4. PBMC Celltypist Labelling

# %% [markdown]
# ### Loops through all 3 loops and collects celltypist predictions -- Models downloaded above

# %%
def celltypist_processing_pbmc(args):
    # Read data
    data_file, model_base_path, output_base_path = args
    pbmc = read_h5_anndata(data_file)

    # Processing testing data
    sample_id = pbmc.obs['pbmc_sample_id'].unique().tolist()[0]
    pbmc.obs.index = pbmc.obs['barcodes']

    # Normalization and log transformation
    sc.pp.normalize_total(pbmc, target_sum=1e4)
    sc.pp.log1p(pbmc)

    # Annotations
    models = {'L1': '18', 'L2': '19', 'L3': '19'}
    predictions = {}

    for level, date in models.items():
        model_file = f'{model_base_path}/ref_pbmc_clean_celltypist_model_flex-features_AIFI_{level}_2024-04-{date}.pkl'
        predictions[level] = celltypist.annotate(pbmc, model=model_file)

        # Write out labels
        predictions[level].predicted_labels.reset_index().to_csv(
            f'{output_base_path}/{sample_id}_{level.lower()}_aifi_predicted_labels.csv')
        # Write out probability matrix
        predictions[level].probability_matrix.reset_index().to_parquet(
            f'{output_base_path}/{sample_id}_{level.lower()}_aifi_probability_matrix.parquet')
        # Write out decision matrix
        predictions[level].decision_matrix.reset_index().to_parquet(
            f'{output_base_path}/{sample_id}_{level.lower()}_aifi_decision_matrix.parquet')


# %% jupyter={"outputs_hidden": true}
args_list = [(file, model_base_path, output_base_path)
             for file in metadata['manual.file_paths']]

with Pool(processes=30) as pool:
    pool.map(celltypist_processing_pbmc, args_list)
