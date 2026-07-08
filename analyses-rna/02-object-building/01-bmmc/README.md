## 01-bmmc-object
#### Code by [Aishwarya Chander](mailto:aishwarya.chander@alleninstitute.org)

Transforms raw BMMC data into a clean, fully annotated single-cell object. Run scripts in numerical order (`01`→`19`); each builds on the previous. Some steps use R.

### Files

- `01-all-build_raw_object.py` — build initial raw BMMC object
- `02-all-first_pass_process.py` — first-pass processing and QC
- `03-all-subject_subsets.py` — subject-specific subsets for plasma-cell analyses
- `04-all-extract_plasma.py` — separate plasma cells via plasma-marker clusters
- `05-plasma-clean_up.ipynb` — refine plasma annotations, remove leaked non-plasma cells
- `06-nonplasma-process.py` — process non-plasma BMMCs
- `07-nonplasma-subsets.py` — L1 cell-type subsets for annotation
- `08-nonplasma-cluster_annotation.ipynb` — annotate clusters with marker genes
- `09-nonplasma-assess_outliers.ipynb` — identify outlier cells
- `10-nonplasma-refine_outliers.ipynb` — separate true outliers from high-expression tumor-associated / progenitor cells
- `11-nonplasma-l1_for_labelling.py` — prepare L1 labels across cell types
- `12-nonplasma-refine_junk_labels.ipynb` — clean low-quality first-pass labels
- `13-all-curate_labels.ipynb` — assemble labels into a single ontology
- `14-all-add_labels.ipynb` — add finalized labels to the raw object
- `15-all-process_objects.ipynb` — final processing of the annotated object
- `16-label-corrections.ipynb` — review and identify label corrections
- `17-all-corrected_curate_labels.ipynb` — re-curate ontology with corrections
- `18-all-corrected_add_labels.ipynb` — add corrected labels to the object
- `19-all-corrected_process_objects.ipynb` — final reprocessing of the corrected object

### Environment

`ac-scanpy-envt` / `ndmm-scrna` (see `../../00-utilities/environment/`); some steps use R.
