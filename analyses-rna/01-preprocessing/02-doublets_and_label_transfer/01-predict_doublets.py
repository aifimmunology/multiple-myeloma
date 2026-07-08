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
# # 2025-07-10: Doublet detection
# ### By [Aishwarya Chander](aishwarya.chander@alleninstitute.org), High Resolution Translational Immunology, Allen Institute for Immunology
#
# **Main aim**: Run scrublet for doublet detection across all samples.

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## 1 Imports and functions

# %% editable=true slideshow={"slide_type": ""}
import hisepy as hp
import os
import pandas as pd
import h5py
import anndata
import scanpy as sc
import scipy.sparse as scs
from concurrent.futures import ThreadPoolExecutor

## Object builders
import sys
sys.path.insert(0, '../../00-utilities/functions/python/')

from h5_data_reader import read_h5_anndata

# %%
file_path_csvs = '../../../data/rna/metadata/'

# %% [markdown]
# ## 2. Get metadata
#
# Pull metadata from HISE to grab all relevant file paths for label prediction. These already have relevant H5AD file paths in them relative to this notebooks path. 

# %% editable=true slideshow={"slide_type": ""}
## Download metadata from HISE using UUIDs
metadata = pd.read_csv('../../../data/rna/metadata/scrna_metadata.csv')


# %% [markdown]
# ## 3. Predict doublets using Scrublet
#
# Using the package `Scrublet` we'll identify and label all the doublets in our dataset. We'll then export it as a .parquet file and call it back while building our batch objects. We'll parallelize these steps across our 235 samples. 

# %%
def scrublet_processing(file_name):
    result = read_h5_anndata(file_name)
    sc.external.pp.scrublet(result)
    return result.obs[['barcodes','predicted_doublet','doublet_score']]


# %% jupyter={"outputs_hidden": true}
# max_workers=30; adjust to match available cores (check with `nproc`)
results = []
with ThreadPoolExecutor(max_workers=30) as executor:  
    for result in executor.map(scrublet_processing, metadata['manual.file_paths']):
        results.append(result)

final_result = pd.concat(results, ignore_index=True)
final_result.to_parquet(file_path_csvs+'all_doublet_scores.parquet')
