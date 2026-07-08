## 02-pbmc-object
#### Code by [Aishwarya Chander](mailto:aishwarya.chander@alleninstitute.org)

Transforms raw PBMC data into clean, labeled single-cell objects (including L1/L2 subsets). Run scripts in numerical order (`01`→`12`); the `04p1`–`04p3` steps are QC-reporting branches off `03` and can be run together.

> ⚠️ 4 samples were added mid-pipeline. QC thresholds, doublet calls, batch integration, and label propagation were re-run so the added samples are handled consistently with the original cohort.

### Files

- `01-all-build-raw_object.ipynb` — build initial raw PBMC object from `.h5` matrices
- `02-all-object_generation.ipynb` — generate the working AnnData object
- `03-all-qc_filtering.ipynb` — QC filtering (gene/count/mito thresholds, doublet removal)
- `04p1-all-qc-doublet_stats.ipynb` — doublet-detection statistics
- `04p2-all-qc-plots-cellfreq.ipynb` — post-filtering cell-frequency QC plots
- `04p3-all-qc-plots-batchqc.ipynb` — batch QC / integration plots
- `05-subset-object_generation.ipynb` — L1 cell-type subset objects
- `06-subset-umaps.ipynb` — per-subset UMAPs
- `07-subset-deepclean.ipynb` — deep-clean subsets (remove leaked/low-quality cells)
- `08-mono_dc-process.py` — process monocyte/DC subset
- `09-mono_dc-correction.ipynb` — correct monocyte/DC annotations
- `10-all-add_labels.ipynb` — add finalized labels to the full PBMC object
- `11-subsample-500k-build_umap.py` — 500k-cell subsampled UMAP for visualization
- `12-all-process-l1-l2-subsets.py` — final processing of L1/L2 subset objects

### Environment

`ac-scanpy-envt` / `ndmm-scrna` (see `../../00-utilities/environment/`).
