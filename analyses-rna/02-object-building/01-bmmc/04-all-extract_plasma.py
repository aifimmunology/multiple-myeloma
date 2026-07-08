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
# # 2025-07-02: Extract BMMC Plasma Cells
# ### By [Aishwarya Chander](aishwarya.chander@alleninstitute.org), High Resolution Translational Immunology, Allen Institute for Immunology
#
# **Main aim**:
# From the subbject specific objects created in the previous notebook, I will now use the label transfer labels to identify approximate neighborhoods of where the Plasma Cells live based on their leiden clusters. If a cluster has any labelled plasma cells, it gets called a plasma cluster in this round. This creates a messy object that mostly exlcudes all non-plasma cells as identified by `CellTypist`.

# %% [markdown]
# ## 1. Imports

# %% editable=true slideshow={"slide_type": ""}
import scanpy as sc
import scanpy.external as sce
import anndata as ad
import pandas as pd

import os
import io
import glob

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec

sc.settings.n_jobs = 30
sc.settings.verbosity = 0

# Set global plotting params
plt.rcParams['figure.dpi'] = 70
plt.rcParams['figure.figsize'] = (5, 5)

# %%
plasma_dir = '../../../data/rna/plasma/'
if not os.path.exists(plasma_dir):
    os.makedirs(plasma_dir)

# %%
ndmm_bmmc_subject_dir = '../../../data/rna/bmmc-subjects/'
all_files = glob.glob(os.path.join(ndmm_bmmc_subject_dir, '*'))
file_paths = [
    f for f in all_files if 'ndmm' not in os.path.basename(f).lower()]

# %% [markdown]
# ## 2. Get data

# %% editable=true slideshow={"slide_type": ""}
anndata_objects = [sc.read_h5ad(file_path) for file_path in file_paths]
len(anndata_objects)

# %% [markdown]
# ## 3. Visualize plasma cell specific markers

# %% editable=true slideshow={"slide_type": ""}
# Initialize list to store plasma cells from each subject
mapped_plasmas = []

# Process each subject's data
for adata in anndata_objects:
    subject = adata.obs['subject.subjectGuid'].unique()[0]

    # Identify cells that are plasma cells from both annotation sources
    adata.obs['plasma'] = ((adata.obs['healthy_l3'] == 'b_plasma') & (
        adata.obs['ext_l3'] == 'Plasma Cell'))
    adata.obs['mapped_plasma'] = False

    # Mark all cells in clusters that contain any plasma cells
    adata.obs.loc[
        adata.obs['leiden'].isin(
            adata.obs[adata.obs['plasma']]['leiden'].unique().tolist()
        ),
        'mapped_plasma'
    ] = True

    # Extract the mapped plasma cells and add to collection
    mapped_plasmas.append(adata[adata.obs['mapped_plasma']].copy())

# Combine and save all plasma cells from all subjects into one dataset
all_plasma = ad.concat(mapped_plasmas, label="subject", index_unique="-")
all_plasma.write(plasma_dir+'messy_raw_plasma.h5ad')

# %%
# Define plasma cell marker genes
plasma_markers = ['CD38', 'TNFRSF17', 'MS4A1', 'CD19', 'JCHAIN', 'TP53', 'MYC', 'KRAS', 'SDC1', 'PRDM1',
                  'IGHG1', 'IGHG2', 'IGHG3', 'IGHG4',
                  'IGHA1', 'IGHA2',
                  'IGHD', 'IGHM', 'IGHE',
                  'IGLC1', 'IGLC2', 'IGLC3', 'IGLC6', 'IGLC7',
                  'IGKC']

# Create PDF file to save all plots for scientist input
if not os.path.exists('results'):
    os.makedirs('results')
    
with PdfPages("../../../data/rna/bmmc/results/plasma_subject_plots.pdf") as pdf:
    # Process each subject's data
    for adata in anndata_objects:
        # Prepare data for analysis
        adata = adata.raw.to_adata()
        sc.pp.normalize_total(adata)  # Normalize to 10k reads per cell
        sc.pp.log1p(adata)  # Log transform
        sc.pp.scale(adata, max_value=10, zero_center=False)  # Scale data

        # Get subject identifier
        subject = adata.obs['subject.subjectGuid'].unique()[0]

        # Create log message with cell counts
        log_buffer = io.StringIO()
        log_buffer.write(f"Cells: {adata.shape[0]} | ")

        # Identify plasma cells from different annotation sources
        adata.obs['ext_plasma'] = adata.obs['ext_l3'] == 'Plasma Cell'
        adata.obs['healthy_bmmc_plasma'] = adata.obs['healthy_l3'] == 'b_plasma'
        ext_count = adata.obs['ext_plasma'].sum()
        healthy_count = adata.obs['healthy_bmmc_plasma'].sum()
        log_buffer.write(
            f"Ext plasma: {ext_count} | Healthy plasma: {healthy_count}")

        # Find cells that are plasma cells in both annotations
        adata.obs['plasma'] = ((adata.obs['healthy_l3'] == 'b_plasma') & (
            adata.obs['ext_l3'] == 'Plasma Cell'))

        # Map all cells in clusters that contain plasma cells
        adata.obs['mapped_plasma'] = False
        adata.obs.loc[adata.obs['leiden'].isin(adata.obs[adata.obs['plasma']]['leiden'].unique().to_list()), 'mapped_plasma'] = (
            True
        )

        # Set up figure layout
        fig = plt.figure(figsize=(18, 16), dpi=70)
        gs = gridspec.GridSpec(
            8, 4, figure=fig,
            height_ratios=[0.05, 0.3, 1.2, 1.15, 1.045, 0.005, 1.4, 1.6],
            hspace=0.25, wspace=0.35
        )

        # Add title
        ax_title = fig.add_subplot(gs[0, :])
        ax_title.axis('off')
        ax_title.text(0.5, 0.5, f"Subject: {subject}", fontsize=22, weight='bold',
                      ha='center', va='center')

        # Add subtitle with cell counts
        ax_subtitle = fig.add_subplot(gs[1, :])
        ax_subtitle.axis('off')
        ax_subtitle.text(0.5, 0.5, log_buffer.getvalue(), fontsize=13, family='monospace',
                         ha='center', va='center')

        # Create UMAP plots for different features
        umap_features = [
            'ext_plasma', 'healthy_bmmc_plasma', 'mapped_plasma', 'leiden',
            'CD38', 'TNFRSF17', 'MS4A1', 'JCHAIN', 'TP53',
            'MYC', 'KRAS', 'SDC1'
        ]
        for i, feat in enumerate(umap_features):
            row = (i // 4) + 2  # Calculate row position
            col = i % 4  # Calculate column position
            ax = fig.add_subplot(gs[row, col])
            sc.pl.umap(adata, color=feat, ax=ax, show=False,
                       title=feat, size=4, frameon=False)
            ax.set_title(feat, fontsize=11, weight='bold')

        # Create dotplot showing marker expression across clusters
        ax_dot = fig.add_subplot(gs[6:8, 0:4])
        sc.pl.dotplot(adata, var_names=plasma_markers, groupby='leiden',
                      ax=ax_dot, show=False, rasterized=True)
        fig.text(
            0.5, ax_dot.get_position().y1 - 0.05,  # Position title above dotplot
            'Plasma Markers Expression by Leiden Cluster',
            ha='center', va='bottom',
            fontsize=14, weight='bold'
        )
        ax_dot.tick_params(axis='x', labelrotation=90, labelsize=9)
        ax_dot.tick_params(axis='y', labelsize=9)

        # Adjust layout and save to PDF
        plt.subplots_adjust(top=0.97, bottom=0.03)
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)
