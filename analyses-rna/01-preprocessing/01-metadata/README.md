## 01-metadata
#### Code by [Aishwarya Chander](mailto:aishwarya.chander@alleninstitute.org)

Sets up base sample metadata for scRNA-seq across BMMC and PBMC and ensures all `.h5` files are downloaded for downstream analysis.

### Files

- `01-make_metadata.ipynb` — builds and formats metadata tables, merging HISE file descriptors with a manually curated clinical metadata sheet (treatment, response, cytogenetics, lab values), deduplicated by specimen ID
- `02-download_files.py` — downloads `.h5` count matrices from HISE and appends local file paths to the metadata CSV

### Output

`scrna_metadata.csv` — sample-level clinical annotations and local file paths.
