# 00-utilities/

Shared functions, reference files, and conda environments for the scRNA-seq pipeline.

## environment/

| File | Purpose |
|---|---|
| `ac-scanpy-envt.yml` | Core Python (scanpy, scrublet, decoupler) |
| `Zi_R_seurat_v5.yml` | Core R (Seurat v5, DESeq2, fgsea, GSVA) |
| `ac-cytetype-envt.yml` | CellTypist cell-type prediction |
| `ac-nmf-envt.yml` | NMF projection (T cell states) |
| `ac-starcat-envt.yml` | starCAT T cell state identification |
| `ac-envt-setup.txt` | Environment creation/activation instructions |
| `scrna_codes.txt` | Cell-type orderings, color palettes, and coding schemes — **use these in all figures; do not redefine inline** |

## functions/python/

| File | Purpose |
|---|---|
| `process_scrna_data.py` | scanpy pipeline: QC, normalization, HVG, PCA, Harmony, UMAP, Leiden |
| `aggregate_counts.py` | Pseudobulk raw-count aggregation (sample × cell type) for DESeq2 |
| `create_frequency_table.py` | Cell-type frequency tables with CLR transformation |
| `filter_genes.py` | Gene filtering by expression thresholds |
| `h5_data_reader.py` | Reads `.h5` count matrices into AnnData |

```python
import sys; sys.path.append('../../00-utilities/functions/python')
from process_scrna_data import process_adata
```

## functions/r/

`base-deseq2.R` — DESeq2 wrappers for all pseudobulk comparisons.

```r
source('../../00-utilities/functions/r/base-deseq2.R')
```
