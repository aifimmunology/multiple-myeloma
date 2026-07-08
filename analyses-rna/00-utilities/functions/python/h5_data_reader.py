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
#     display_name: ''
#     name: ''
# ---

# %%
"""
Readers for AIFI .h5 scRNA-seq files.

Constructs AnnData objects from .h5 files produced by the HISE pipeline.
Used across doublet detection and CellTypist label prediction steps.
"""

import h5py
import pandas as pd
import scipy.sparse as scs
import anndata


def read_expression_matrix(h5_file):
    """Reads the sparse expression matrix from an open h5 file handle."""
    matrix = scs.csc_matrix(
        (h5_file['matrix']['data'][:],
         h5_file['matrix']['indices'][:],
         h5_file['matrix']['indptr'][:]),
        shape=tuple(h5_file['matrix']['shape'][:])
    )
    return matrix


def read_metadata(h5_file):
    """Reads per-cell metadata (barcodes + observation columns) from an open h5 file handle."""
    barcodes = [b.decode('UTF-8') for b in h5_file['matrix']['barcodes'][:]]
    metadata_df = pd.DataFrame({'barcodes': barcodes})

    for column in h5_file['matrix']['observations'].keys():
        values = h5_file['matrix']['observations'][column][:]
        if isinstance(values[0], (bytes, bytearray)):
            values = [v.decode('UTF-8') for v in values]
        metadata_df[column] = values

    return metadata_df


def read_h5_anndata(h5_file_path):
    """
    Constructs an AnnData object from an AIFI .h5 file.

    Parameters
    ----------
    h5_file_path : str
        Path to the .h5 count matrix file.

    Returns
    -------
    anndata.AnnData
        AnnData with cells x genes, barcodes in obs, gene names in var.
    """
    with h5py.File(h5_file_path, mode='r') as h5_file:
        expression_matrix = read_expression_matrix(h5_file)
        gene_names = [g.decode('UTF-8') for g in h5_file['matrix']['features']['name'][:]]
        metadata_df = read_metadata(h5_file)

        adata = anndata.AnnData(expression_matrix.T, obs=metadata_df)
        adata.var_names = gene_names
        adata.var_names_make_unique()

    return adata
