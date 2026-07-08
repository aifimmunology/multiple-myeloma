# Standard library imports
import os
import fnmatch
from IPython.display import display, Markdown
import difflib
from datetime import date
import re
import time
import gc

# Third party imports
import pandas as pd 
import numpy as np 
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import anndata as ad

# # Specific module imports
# import hisepy

# # Ignore warnings
# import warnings
# warnings.filterwarnings('ignore')
# warnings.filterwarnings('ignore', category=FutureWarning, module='seaborn._oldcore')

# Scanpy settings
sc.settings.n_jobs = 30
sc.settings.verbosity = 0
sc.settings.set_figure_params(
    dpi=80,
    facecolor='white',
    frameon=False,
    fontsize=10, 
    dpi_save=400,
    figsize=(6,6), 
    format='png'
)

def get_ram_usage():
    # Getting all memory using os.popen()
    total_memory, used_memory, free_memory = map(int, os.popen('free -t -m').readlines()[-1].split()[1:])
    # Memory usage
    print("RAM memory % used:", round((used_memory / total_memory) * 100, 2))
    print(f"total_memory = {total_memory}, used_memory = {used_memory}, free_memory = {free_memory}")

def get_hise_uuids(project_store_name, destination_tag):
    all_uuids = hisepy.list_files_in_project_store(str(project_store_name))
    file_identifier = [str(destination_tag)]
    file_uuids_df = all_uuids[all_uuids['name'].apply(lambda x: any(i in x for i in file_identifier))][['id', 'name']]
    h5ad_uuids = file_uuids_df['id'].to_list()
    return file_uuids_df, h5ad_uuids

def find_file_path_with_name(start_dir, target_file):
    for root, _, files in os.walk(start_dir):
        if target_file in files:
            return os.path.join(root, target_file)
    return None

def find_file_paths_with_extension(directory, extension):
    matches = []
    for root, dirnames, filenames in os.walk(directory):
        for filename in fnmatch.filter(filenames, '*.' + extension):
            matches.append(os.path.join(root, filename))
    return matches

def count_files_in_directory(directory_path):
    return len([name for name in os.listdir(directory_path) if os.path.isfile(os.path.join(directory_path, name))])

def display_hise_params():
    code = '''
    ## Create a list to add all the files to download from HISE
    download_files = []
    
    ## Create a list to add all the files to upload to HISE
    upload_files = []
    
    ## Define destination ID
    import random, string
    destination_tag = ''.join(random.choices(string.ascii_letters + string.digits, k=10))
    
    ## Study space UUID
    colab_space_name = hp.get_study_spaces()[6]['name']
    study_space_uuid = hp.get_study_spaces()[6]['id']
    
    ## Upload title
    from datetime import date
    title = f'{str(date.today())} Add title here [AC]'
    '''
    display(Markdown(f'```python\n{code}\n```'))

def check_for_genes(adata, genes):
    var_names = adata.raw.var_names.str.upper()
    
    for gene in genes:
        gene_upper = gene.upper()
        if gene_upper in var_names:
            print(f"{gene}: Exists, {gene_upper}")
        else:
            closest_match = difflib.get_close_matches(gene_upper, var_names, n=1)
            if closest_match:
                print(f"{gene}: Maybe try, {closest_match[0]}")
            else:
                print(f"{gene}: No close match found")

def filter_markers(markers_dict, valid_markers):
    # Filter the dictionary
    filtered_dict = {k: [marker for marker in v if marker in valid_markers] for k, v in markers_dict.items()}
    # Remove keys with empty lists
    filtered_dict = {k: v for k, v in filtered_dict.items() if v}
    return filtered_dict

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

def assign_labels_by_leiden(adata, label_key, leiden_key="leiden"):
    # Create a DataFrame with cluster assignments and labels
    df = pd.DataFrame({
        'cluster': adata.obs[leiden_key],
        'label': adata.obs[label_key]
    })
    
    # Find majority label for each cluster
    majority_labels = (df[df['label'].notna()]
                      .groupby('cluster')['label']
                      .agg(lambda x: x.mode().iloc[0] if len(x.mode()) > 0 else None)
                      .to_dict())
    
    # Create new column with propagated labels
    adata.obs[label_key + '_leiden'] = adata.obs[label_key].copy()
    mask = adata.obs[label_key].isna()
    adata.obs.loc[mask, label_key + '_leiden'] = adata.obs.loc[mask, leiden_key].map(majority_labels)


def rank_genes_and_plot(adata, 
                        groupby, 
                        method='t-test', 
                        corr_method='benjamini-hochberg', 
                        n_genes=4, 
                        cmap='bwr', 
                        min_logfoldchange=2):
    # Split the groupby string by '.' and take the last segment
    key_segment = groupby.split('.')[-1]
    
    # Construct the key and key_added
    key_added = f'{key_segment}-{method}'
    
    # Perform rank genes groups analysis
    sc.tl.rank_genes_groups(
        adata=adata,
        groupby=groupby,
        method=method,
        corr_method=corr_method,
        key_added=key_added
    )
    
    # Compute the dendrogram
    sc.tl.dendrogram(adata, groupby=groupby)
    
    # Plot the rank genes groups dotplot
    sc.pl.rank_genes_groups_dotplot(
        adata, 
        n_genes=n_genes, 
        key=key_added, 
        groupby=groupby,     
        values_to_plot='logfoldchanges',
        cmap=cmap,
        min_logfoldchange=min_logfoldchange
    )

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


def get_labels_and_scores(pbmc_sample_id, data_sources, label_dir, levels=3):
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
        label_dfs = [pd.read_csv(f'{label_dir}/{pbmc_sample_id}_l{i}_{data_source}_predicted_labels.csv')[['barcodes', 'predicted_labels']] for i in range(1, levels+1)]
        
        merged_label_df = reduce(lambda left,right: pd.merge(left,right,on='barcodes', how='left'), label_dfs)
        merged_label_df.columns = ['barcodes'] + [f'{data_source}_l{i}' for i in range(1, levels+1)]
        
        return merged_label_df

    merged_labels = [load_and_merge(data_source) for data_source in data_sources]
    
    return reduce(lambda left,right: pd.merge(left,right,on='barcodes', how='left'), merged_labels)


def add_metadata(adata, metadata):
    """Adds metadata to an AnnData object."""
    pbmc_sample_id = adata.obs['pbmc_sample_id'][0]
    df = get_labels_and_scores(pbmc_sample_id, ['ext', 'healthy', 'aifi'])

    # define where to find the doublet scores
    doublet_scores = pd.read_parquet('../output/ndmm-bmmc-parquet/2025-01-24_all_bmmc_doublet_scores.parquet')
    df = df.merge(doublet_scores, on='barcodes', how='left')
    adata.obs = adata.obs.merge(df, on='barcodes', how='left')

    # what metadata you want your adata to have
    add_meta_cols = [
        'pbmc_sample_id', 
        'sample.sampleKitGuid', 
        'sample.visitDetails',
        'sample.drawDate', 
        'sample.daysSinceFirstVisit',
        'sample.diseaseStatesRecordedAtVisit', 
        'subject.biologicalSex', 
        'subject.birthYear', 
        'subject.ethnicity',
        'subject.partnerCode', 
        'subject.race', 
        'subject.subjectGuid', 
        'cohort.cohortGuid',
        'manual.time_stamp',
        'manual.tissue', 
        'manual.disease_condition', 
        'manual.response',
        'manual.response_type', 
        'manual.extracted_name', 
        'manual.file_paths',
        'manual.batch_id',
        'manual.cmv.ab_screen_index_value',
        'manual.cmv.ab_screen_result', 
        'manual.cmv.ab_visit', 
        'manual.igg.iga',
        'manual.igg.igg', 
        'manual.igg.igm', 
        'manual.igg.immunofixation',
        'manual.igg.k_flc', 
        'manual.igg.kl_flc_ratio', 
        'manual.igg.l_flc',
        'manual.igg.mono_1_id', 
        'manual.igg.mono_1_quant',
        'manual.lip.cholesterol_hdl_ratio', 
        'manual.lip.cholesterol_hdl',
        'manual.lip.cholesterol_ldl', 
        'manual.lip.cholesterol_non_hdl',
        'manual.lip.cholesterol_total', 
        'manual.lip.triglycerides',
        'manual.RISS.stage', 
        'manual.cyto.17p_loss', 
        'manual.cyto.1q_gain',
        'manual.cyto.1p_loss', 
        'manual.cyto.high_risk', 
        'manual.cyto.karyotype',
        'manual.cyto.t11_14',
        'manual.cyto.t14_16', 
        'manual.cyto.t14_20',
        'manual.cyto.t4_14', 
        'manual.cyto.karyotype.1']
    
    add_meta = metadata[add_meta_cols]
    adata.obs = adata.obs.merge(add_meta, on='pbmc_sample_id', how='left')
    adata.obs = adata.obs.drop('barcodes', axis = 1)
    return adata

def build_adata(h5_file, metadata, output_dir, save=False):
    """Builds an AnnData object from an HDF5 file."""
    h5_con = h5py.File(h5_file, mode='r')
    rna_mat = read_mat(h5_con, 'matrix')
    adt_mat = read_mat(h5_con, 'ADT')
    obs = read_obs(h5_con)
    obs = obs.reset_index(drop=True)
    barcodes = obs['barcodes']
    obs = obs.drop('barcodes', axis=1)
    genes = read_feats(h5_con, 'matrix', 'name')
    adts = read_feats(h5_con, 'ADT', 'id')
    h5_con.close()
    
    adata = sc.AnnData(X=rna_mat.T, obs=obs)
    adata.var_names = genes
    adata.var_names_make_unique()
    adata.obs_names = barcodes
    
    adata = add_metadata(adata, metadata)
    adata.obs_names = barcodes
    
    adt_mat = adt_mat.toarray().T
    adt_df = pd.DataFrame(adt_mat, columns=adts, index=barcodes)
    adata.obsm['adt_counts'] = adt_df

    pbmc_sample_id = adata.obs['pbmc_sample_id'][0]
    
    if save:
        adata.write_h5ad(f'{output_dir}/adata_raw_bmmc_{pbmc_sample_id}.h5ad')
    
    return adata
