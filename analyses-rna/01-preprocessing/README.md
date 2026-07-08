## 01-preprocessing

Metadata preparation, doublet detection, and cell-type prediction for the scRNA-seq cohort.

### 01-metadata/

| File | Description |
|---|---|
| `01-make_metadata.ipynb` | Builds sample metadata from HISE descriptors, merges curated clinical annotations (treatment, response, cytogenetics, labs), deduplicates by specimen ID |
| `02-download_files.py` | Downloads `.h5` count matrices from HISE and appends local paths to the metadata CSV |

Output: `scrna_metadata.csv` (sample-level clinical + file-path info).

### 02-doublets_and_label_transfer/

| File | Description |
|---|---|
| `01-predict_doublets.py` | Per-sample doublet detection (Scrublet) |
| `02-predict_labels.py` | Cell-type prediction (CellTypist) against 3 references |

Output: per-sample doublet scores and cell-type predictions.

### Run order

`01-metadata/01-make_metadata.ipynb → 01-metadata/02-download_files.py → 02-doublets_and_label_transfer/01-predict_doublets.py → 02-predict_labels.py`

### Environments

`ac-scanpy-envt.yml` (Python) throughout; `Zi_R_seurat_v5.yml` for the R metadata step.
