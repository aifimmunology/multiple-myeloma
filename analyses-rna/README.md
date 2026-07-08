# analyses-rna/

Single-cell RNA-seq analysis of BMMC and PBMC from the NDMM cohort: raw count matrices → cell-type annotation → subset analyses → pseudobulk differential expression.

## Structure

```
00-utilities/       # Shared functions, environments, reference files
01-preprocessing/   # Metadata, doublet detection, cell-type prediction
02-object-building/ # AnnData construction (BMMC, PBMC)
03-subset-analyses/ # Plasma cell, T cell, monocyte/DC deep dives
04-pseudobulk/      # Pseudobulk aggregation, DESeq2, FGSEA, GSVA
```

`04-pseudobulk/` contains the primary `l3-pseudobulk/` pipeline (frequencies → DESeq2 → FGSEA) plus peer-review additions: `revision-gsva/` (per-sample GSVA) and Wald-statistic-ranked FGSEA (`revision-walds-fgsea-no-filter.R`, `walds_eval.ipynb`).

## Environments

| Environment | Config | Use |
|---|---|---|
| `ndmm-scrna` | `00-utilities/environment/ac-scanpy-envt.yml` | Core Python (scanpy, scrublet, decoupler) |
| `r-seurat` | `00-utilities/environment/Zi_R_seurat_v5.yml` | Core R (Seurat v5, DESeq2, fgsea, GSVA) |
| `ac-cytetype-envt` | `00-utilities/environment/ac-cytetype-envt.yml` | CellTypist prediction |
| `ac-nmf-envt` | `00-utilities/environment/ac-nmf-envt.yml` | NMF projection (T cell states) |
| `ac-starcat-envt` | `00-utilities/environment/ac-starcat-envt.yml` | starCAT T cell states |

Setup: `00-utilities/environment/ac-envt-setup.txt`.

## Run order

`01-preprocessing/ → 02-object-building/ → 03-subset-analyses/ → 04-pseudobulk/`. Scripts within each directory are numbered; see subdirectory READMEs.

## Data

Objects, labels, CellTypist predictions, and pseudobulk/DE tables under `data/rna/`; gene sets under `data/gmt/`.

## Contributors

- InferCNV analyses by [Palak Genge](mailto:palak.genge@alleninstitute.org)
- PBMC annotation cleanup by [Mansi Singh](mailto:mansi.singh@alleninstitute.org)
- All other analyses by [Aishwarya Chander](mailto:aishwarya.chander@alleninstitute.org)
- Thanks to [Lucas Graybuck](mailto:lucas.graybuck@alleninstitute.org), [Ziyuan He](mailto:ziyuan.he@alleninstitute.org), and [Qiuyu Gong](mailto:qiuyu.gong@alleninstitute.org) for scRNA-seq expertise.
