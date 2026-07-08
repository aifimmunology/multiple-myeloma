## 04-pseudobulk
#### Code by [Aishwarya Chander](mailto:aishwarya.chander@alleninstitute.org)

L3 cell-type pseudobulk analysis: cell-type frequencies, DESeq2 differential expression (longitudinal and vs. healthy), FGSEA pathway enrichment, and — added during peer review — per-sample GSVA enrichment. BMMC and PBMC are run as parallel `bmmc`/`pbmc` script pairs.

### l3-pseudobulk/

Primary pipeline; run in numerical order.

| Stage | Files | Description |
|---|---|---|
| Labels | `01-add-aifi_plot_l3-labels.ipynb` | Attach L3 labels / metadata to the final objects |
| Frequencies | `02-create-frequencies.ipynb`, `03-run-timepoint_freq_analysis.ipynb`, `04`/`05-plot-*-timepoint-boxplots.ipynb` | CLR cell-type frequencies, timepoint comparison stats, boxplots |
| Pseudobulk prep | `06-bmmc-deseq2-prep.ipynb`, `07-pbmc-deseq2-prep.ipynb` | Sum raw counts per sampleKit × L3 cell type; expressed-gene filtering; sample-kit covariate tables |
| DESeq2 | `08`/`09-*-time_points-deseq2.ipynb`, `10`/`11-*-timepoints_vs_healthy-deseq2.ipynb` | Differential expression (longitudinal and vs. healthy) |
| DE plots | `12-plot-deseq2-bars.R`, `13-plots-deseq2-volcano_fixed.ipynb` | DE-gene count bars and volcano plots |
| Pathways | `14-get_gmt_pathways.R`, `15-all-fgsea-no-filter.R`, `16-plot-fgsea-no-filter.ipynb` | Download gene sets, run FGSEA, plot dotplots |
| Revision | `revision-walds-fgsea-no-filter.R`, `walds_eval.ipynb` | Wald-statistic-ranked FGSEA and ranking-metric evaluation |

Gene-set databases live in `gmt/`.

### revision-gsva/  (peer-review addition)

Per-sample GSVA enrichment: `01-create_filtered_gene_lists` → `02`/`03-*_run_gsva` (bmmc/pbmc) → `04`/`05-*_significance` (pbmc/bmmc) → `06-relax_pvals` → `07-plot_violins`.

### Methods

- **Pseudobulk:** summed raw counts per sample (`sampleKitGuid`) × L3 cell type; genes filtered to the expressed set per contrast before DESeq2.
- **DESeq2 designs:**
  - Longitudinal (timepoint vs. timepoint): `~ subject.subjectGuid + label.visitDetails` (within-subject / paired).
  - Cross-sectional (vs. healthy): `~ label.visitDetails + subject.biologicalSex + subject.age + subject.cmv`.
- **FGSEA:** run per cell type on the DESeq2 gene rankings; a peer-review variant ranks genes by the DESeq2 Wald statistic (`walds_eval.ipynb` compares ranking metrics). Enriched pathways reported at adjusted P < 0.05.
- **Frequencies:** CLR-transformed cell-type frequencies, compared across timepoints in `03`.

### Environment

Python (`ac-scanpy-envt` / `ndmm-scrna`) for prep, frequencies, and GSVA; R (`r-seurat`: DESeq2, fgsea, GSVA) for differential expression and enrichment.

### Data

Inputs: final labelled objects from `02-object-building/` (`data/rna/`). Outputs: pseudobulk matrices and DESeq2 / FGSEA / GSVA result tables under `data/rna/pseudobulk/` and `data/rna/gsva/`; gene sets under `data/gmt/`.
