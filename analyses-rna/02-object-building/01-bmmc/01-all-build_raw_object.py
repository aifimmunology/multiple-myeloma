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
# # 2025-07-01: Longitudinal MM Data build
# ### By [Aishwarya Chander](aishwarya.chander@alleninstitute.org), High Resolution Translational Immunology, Allen Institute for Immunology
#
# **Main aim**: 
# In this notebook, I take my previously cleaned up metadata, predicted doublets and CellTypist lables and add them onto an object that contains all BMMC samples. This will act as my base object for further analyses. I also save all objects per sample. 

# %% [markdown]
# ## 1. Imports

# %%
import os
import warnings
from functools import reduce

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import scanpy as sc
import anndata
import h5py
import scipy.sparse as scs
import concurrent.futures
from tqdm import tqdm

import hisepy as hp

# Settings
warnings.filterwarnings("ignore")
sc.settings.n_jobs = 60
plt.style.use('default')

print("Current working directory:", os.getcwd())

# %% [markdown]
# #### Necessary Paths

# %% editable=true slideshow={"slide_type": ""}
celltypist_labels_dir = '../../../data/rna/celltypist/labels/'
file_path_csvs = '../../../data/rna/metadata/'

bmmc_subject_dir = '../../../data/rna/bmmc-subjects/'
if not os.path.exists(bmmc_subject_dir):
    os.makedirs(bmmc_subject_dir)


# %% [markdown] editable=true slideshow={"slide_type": ""}
# ### 1.2. Custom Functions
# > It ain't much but its honest work. 

# %% editable=true slideshow={"slide_type": ""}
def check_list_and_add_files(given_list, files):
    """
    Adds files to the given list if they are not already present.

    Parameters:
    given_list (list): The list to which files will be added.
    files (list or str): The files to be added. Can be a single file (str) or a list of files (list).
    """
    if isinstance(files, str):
        files = [files]
    
    for file in files:
        if file not in given_list:
            given_list.append(file)
            
def read_mat(h5_con, mat_name):
    """Reads a matrix from an HDF5 file."""
    mat = scs.csc_matrix(
        (h5_con[mat_name]['data'][:], 
         h5_con[mat_name]['indices'][:], 
         h5_con[mat_name]['indptr'][:]), 
        shape = tuple(h5_con[mat_name]['shape'][:])
    )
    return mat

def read_feats(h5_con, mat_name, name_col):
    """Reads features from an HDF5 file."""
    feats = h5_con[mat_name]['features'][name_col][:]
    feats = [x.decode('UTF-8') for x in feats]
    return feats

def read_obs(h5con):
    """Reads observations from an HDF5 file."""
    bc = h5con['matrix']['barcodes'][:]
    bc = [x.decode('UTF-8') for x in bc]
    obs_df = pd.DataFrame({ 'barcodes' : bc })
    obs_columns = h5con['matrix']['observations'].keys()
    for col in obs_columns:
        values = h5con['matrix']['observations'][col][:]
        if(isinstance(values[0], (bytes, bytearray))):
            values = [x.decode('UTF-8') for x in values]
        obs_df[col] = values
    return obs_df

def read_scrublet(scrublet_file):
    """Reads a scrublet file."""
    scrub = pd.read_csv(scrublet_file, index_col = 0)
    scrub.index = scrub['barcodes']
    scrub = scrub.drop('barcodes', axis = 1)
    return scrub

def read_labels(labels_file):
    """Reads a labels file."""
    labels = pd.read_csv(labels_file)
    labels.index = labels['barcodes']
    labels = labels.drop('barcodes', axis = 1)
    return labels


# %% editable=true slideshow={"slide_type": ""}
## Checkpoint -- Make sure path is correct
celltypist_labels_dir = '../../../data/rna/celltypist/labels/'

def get_labels_and_scores(pbmc_sample_id, data_sources, celltypist_labels_dir=celltypist_labels_dir, levels=3):
    """
    For each pbmc_sample_id, extract labels predicted by different methods and combine them into a DataFrame.

    Parameters:
    pbmc_sample_id (str): The sample id to process.
    data_sources (list): List of data sources for different methods.
    levels (int): Number of levels of labels to process. Default is 3.

    Returns:
    DataFrame: DataFrame with the combined labels.
    """
    def load_and_merge(data_source):
        label_dfs = [pd.read_csv(f'{celltypist_labels_dir}{pbmc_sample_id}_l{i}_{data_source}_predicted_labels.csv')[['barcodes', 'predicted_labels']] for i in range(1, levels+1)]
        
        merged_label_df = reduce(lambda left,right: pd.merge(left,right,on='barcodes', how='left'), label_dfs)
        merged_label_df.columns = ['barcodes'] + [f'{data_source}_l{i}' for i in range(1, levels+1)]
        
        return merged_label_df

    merged_labels = [load_and_merge(data_source) for data_source in data_sources]
    
    return reduce(lambda left,right: pd.merge(left,right,on='barcodes', how='left'), merged_labels)


# %%
## Checkpoint -- Make sure path is correct
def add_metadata(adata, metadata, doublet_file_path):
    """Adds metadata to an AnnData object."""
    pbmc_sample_id = adata.obs['pbmc_sample_id'][0]
    df = get_labels_and_scores(pbmc_sample_id, ['ext', 'healthy', 'aifi'])

    # define where to find the doublet scores
    doublet_scores = pd.read_parquet(doublet_file_path)
    df = df.merge(doublet_scores, on='barcodes', how='left')
    adata.obs = adata.obs.merge(df, on='barcodes', how='left')

    # what metadata you want your adata to have
    add_meta_cols = [
        'pbmc_sample_id', 
        'sample.sampleKitGuid', 
        'sample.visitDetails',
        'sample.visitName',
        'sample.drawDate', 
        'sample.daysSinceFirstVisit',
        'sample.diseaseStatesRecordedAtVisit', 
        'subject.biologicalSex', 
        'subject.birthYear', 
        'subject.ethnicity',
        'subject.partnerCode', 
        'subject.race', 
        'subject.subjectGuid', 
        'subject.cmv',
        'specimen.specimenGuid', 
        'cohort.cohortGuid',
        'manual.time_stamp',
        'tissue', 
        'manual.response',
        'manual.response_type', 
        'manual.extracted_name', 
        'manual.batch_id', 
        'manual.category', 
        'manual.treatment_dara', 
        'manual.flu_response']
    
    add_meta = metadata[add_meta_cols]
    adata.obs = adata.obs.merge(add_meta, on='pbmc_sample_id', how='left')
    adata.obs = adata.obs.drop('barcodes', axis = 1)
    return adata


# %%
## Checkpoint -- Make sure path is correct
bmmc_raw_dir = '../../../data/rna/raw-files/'
if not os.path.exists(bmmc_raw_dir):
    os.makedirs(bmmc_raw_dir)

def build_adata(h5_file, metadata, doublet_file_path, bmmc_raw_dir=bmmc_raw_dir, save=False):
    """Builds an AnnData object from an HDF5 file."""
    h5_con = h5py.File(h5_file, mode = 'r')
    rna_mat = read_mat(h5_con, 'matrix')
    adt_mat = read_mat(h5_con, 'ADT')
    obs = read_obs(h5_con)
    obs = obs.reset_index(drop = True)
    barcodes = obs['barcodes']
    obs = obs.drop('barcodes', axis = 1)
    genes = read_feats(h5_con, 'matrix', 'name')
    adts = read_feats(h5_con, 'ADT', 'id')
    h5_con.close()
    adata = sc.AnnData(X = rna_mat.T, obs = obs)
    adata.var_names = genes
    adata.var_names_make_unique()
    adata.obs_names = barcodes
    
    adata = add_metadata(adata, metadata, doublet_file_path)
    adata.obs_names = barcodes
    
    adt_mat = adt_mat.toarray().T
    adt_df = pd.DataFrame(adt_mat, columns = adts, index = barcodes)
    adata.obsm['adt_counts'] = adt_df

    if save:
        pbmc_sample_id = adata.obs['pbmc_sample_id'][0]
        adata.write_h5ad(f'{bmmc_raw_dir}adata_raw_bmmc_{pbmc_sample_id}.h5ad')
    
    return adata


# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## 2. Grab files

# %% [markdown]
# ### 2.1. Get Metadata

# %%
metadata = pd.read_csv('../../../data/rna/metadata/scrna_metadata.csv')
metadata = metadata[(metadata['tissue'] == 'BMMC') & (metadata['manual.treatment_dara'] == 'non_dara')]
metadata = metadata.fillna('None')

file_paths = metadata['manual.file_paths'].to_list()
metadata = metadata.drop(columns=['Unnamed: 0', 'Unnamed: 0.1'], errors='ignore')

# %%
doublet_uuid = ['146ab3d0-fb75-4d22-9817-a87009921a45']
doublet_file_path = hp.cache_files(doublet_uuid)
doublet_file_path = doublet_file_path[0]


# %%
def count_files_in_directory(directory_path):
    return len([name for name in os.listdir(directory_path) if os.path.isfile(os.path.join(directory_path, name))])

file_path_labels = '../../../data/rna/celltypist/labels/'
file_count = count_files_in_directory(file_path_labels)
print(f'There are {file_count} files in the directory.')

# %%
files = list(metadata['file.id'])
len(files)

# %%
file_paths = metadata['manual.file_paths'].to_list()

# %% [markdown] editable=true slideshow={"slide_type": ""}
# ## 3. Build and save master object

# %%
from concurrent.futures import ThreadPoolExecutor, as_completed

# %%
metadata

# %%
adata_list = []

with ThreadPoolExecutor(max_workers=20) as executor:
    future_to_file = {executor.submit(build_adata, file_path, metadata, doublet_file_path, bmmc_raw_dir=bmmc_raw_dir, save=True): file_path for file_path in file_paths}
    for future in tqdm(as_completed(future_to_file), total=len(file_paths)):
        result = future.result()
        if result is not None:
            adata_list.append(result)

# %%
adata = anndata.concat(adata_list)
adata.uns['all_sample_metadata'] = metadata.set_index('pbmc_sample_id').to_dict(orient='index')

# %%
adata.write_h5ad(bmmc_raw_dir+'all-bmmc-raw.h5ad')
adata.obs.to_csv(bmmc_raw_dir+'all-bmmc-metadata.csv')
