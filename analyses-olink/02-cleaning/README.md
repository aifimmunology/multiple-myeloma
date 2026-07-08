# Olink object generation

Builds two Olink objects from raw HISE batches: **BMIF** (NDMM) and **Plasma** (NDMM + selected BR1/BR2 healthy donors). Cleaning is serialized (`01`→`06`) and:

1. Excludes daratumumab-treated subjects.
2. Includes selected BR1/BR2 healthy donors (Plasma object only).
3. Retains only NDMM treatment and flu-series timepoints.

## Dependencies

Uses functions from the internal packages `ndmmFH1` and `olink-qc-report`; where needed, those functions are copied into and sourced from this repository.

## Source data (HISE)

Bone marrow — `01-bmif-pull-bm-samples.ipynb`:
- `Q-11665_Skene_NPX_2023-12-05.xlsx` — `628523de-7712-4668-a93b-787d859be647`

Plasma — `04-pl-pull-plasma-samples.ipynb` (8 bridged batches):
- `00753770-4803-4947-a290-e7469c410067`
- `1a528317-bbf3-45b7-82c9-529577cf0b15`
- `55d5d20a-506a-459e-8315-7d66276fb8f9`
- `59f5f656-085f-4d0b-a240-68a9eb15e68e`
- `65103a13-52de-4f65-89fa-dd03296b4ecb`
- `7c2d668a-b851-4087-ba31-408de4ce889c`
- `b494cce8-1314-4f6e-9666-42f0c6e1c702`
- `748f218f-07c9-4303-a500-97b7621f920c`

Sample metadata is joined from HISE file descriptors (`hise::getFileDescriptors(fileType = "Olink")`).

## Output

`06-pl-bmif-align-final-olink-objects.ipynb` produces `MM_BMIF_Olink_Final.RDS` and `MM_Plasma_Olink_Final.RDS`, the basis for all downstream analyses.
