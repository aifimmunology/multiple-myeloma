import pandas as pd
import numpy as np
import scanpy as sc
import os
import re

# Reference: https://github.com/aifimmunology/IHA-Figure/blob/79ea9445825dc9c82e3bf07bba719a7c1e0bbe6d/helper_function/helper_function_IHA.py


def grouped_obs_sum_raw(adata, group_key, layer=None, gene_symbols=None):
    '''
    Compute the sum of raw counts for each group defined by `group_key` in the given AnnData object.

    Parameters:
    -----------
    adata : AnnData
        Filtered AnnData object containing the single-cell data.
    group_key : str
        Key in `adata.obs` to group the observations by.
    layer : str, optional
        If specified, the layer in `adata` to use for calculations. If None, use `adata.X`.
    gene_symbols : list of str, optional
        List of gene symbols to include in the sum. If None, include all genes.

    Returns:
    --------
    pd.DataFrame
        DataFrame containing the summed observations for each group. Rows correspond to genes and columns to groups.
    '''
    if layer is not None:
        def getX(x): return x.layers[layer]
    else:
        def getX(x): return x.X

    if gene_symbols is not None:
        idx = adata.var_names.isin(gene_symbols)
        new_idx = adata.var_names[idx]
    else:
        new_idx = adata.var_names

    grouped = adata.obs.groupby(group_key, observed=False)
    out = pd.DataFrame(
        np.zeros((len(new_idx), len(grouped)), dtype=np.float64),
        columns=list(grouped.groups.keys()),
        index=new_idx
    )

    for group, idx in grouped.indices.items():
        X = getX(adata[idx])
        out[group] = np.ravel(X.sum(axis=0, dtype=np.float64))

    return out


def grouped_obs_mean(adata, group_key, layer=None, gene_symbols=None):
    '''
    Compute the mean of log-normalized expression for each group defined by `group_key` in the given AnnData object.

    Parameters:
    -----------
    adata : AnnData
        Filtered AnnData object containing the single-cell data.
    group_key : str
        Key in `adata.obs` to group the observations by.
    layer : str, optional
        If specified, the layer in `adata` to use for calculations. If None, use `adata.X`.
    gene_symbols : list of str, optional
        List of gene symbols to include in the mean calculation. If None, include all genes.

    Returns:
    --------
    pd.DataFrame
        DataFrame containing the mean observations for each group. Rows correspond to genes and columns to groups.
    '''
    if layer is not None:
        def getX(x): return x.layers[layer]
    else:
        def getX(x): return x.X

    if gene_symbols is not None:
        idx = adata.var_names.isin(gene_symbols)
        new_idx = adata.var_names[idx]
    else:
        new_idx = adata.var_names

    grouped = adata.obs.groupby(group_key, observed=False)
    out = pd.DataFrame(
        np.zeros((len(new_idx), len(grouped)), dtype=np.float64),
        columns=list(grouped.groups.keys()),
        index=new_idx
    )

    for group, idx in grouped.indices.items():
        X = getX(adata[idx])
        out[group] = np.ravel(X.mean(axis=0, dtype=np.float64))

    return out


def process_sample(sample_id, adata, group_key, results_dir, tissue):

    match = re.search(r'(l\d)', group_key)
    if not match:
        raise ValueError(f"Could not extract level id from {group_key}")
    level_id = match.group(1)

    # Ensure output directories exist
    os.makedirs(
        f'{results_dir}/{tissue}_{level_id}_raw_gexp_celltypes_per_samplekit', exist_ok=True)
    os.makedirs(
        f'{results_dir}/{tissue}_{level_id}_lognorm_gexp_celltypes_per_samplekit', exist_ok=True)

    # Compute raw sum
    raw_sum_df = grouped_obs_sum_raw(adata, group_key)
    raw_sum_df.to_csv(
        f'{results_dir}/{tissue}_{level_id}_raw_gexp_celltypes_per_samplekit/{sample_id}.csv')

    # Normalize and log-transform
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Compute mean
    mean_df = grouped_obs_mean(adata, group_key)
    mean_df.to_csv(
        f'{results_dir}/{tissue}_{level_id}_lognorm_gexp_celltypes_per_samplekit/{sample_id}.csv')


def process_sample_wrapper(args):
    return process_sample(*args)
