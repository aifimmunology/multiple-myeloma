# analyses-olink/

Olink proteomics cleaning and longitudinal analysis for the NDMM cohort, across two sample types: bone marrow interstitial fluid (BMIF) and plasma.

## Structure

```
00-utilities/   # Shared functions and environment files
01-analysis/    # Longitudinal, cross-sectional, and pathway analyses
02-cleaning/    # Data ingestion, batch alignment, object construction
```

## Run order

Run `02-cleaning/` first to build the final Olink objects, then `01-analysis/`.

**02-cleaning/**

| Script | Output |
|---|---|
| `01-bmif-pull-bm-samples.ipynb` | `BM_Olink_Samples.csv` |
| `02-bmif-join-metadata.ipynb` | `Olink_BM_with_hise_descriptors.csv` |
| `03-bmif-clean-dataset.ipynb` | `Data_Olink_BM_cleaned.csv` |
| `04-pl-pull-plasma-samples.ipynb` | combined bridged plasma CSV |
| `05-pl-join-metadata.ipynb` | `Olink_Plasma_FH1_BR1_BR2_with_hise_descriptors.csv` |
| `06-pl-bmif-align-final-olink-objects.ipynb` | `MM_BMIF_Olink_Final.RDS`, `MM_Plasma_Olink_Final.RDS` |

**01-analysis/**

| Script | Description |
|---|---|
| `01-longitudinal_pairwise_analyses.ipynb` | Pairwise comparisons across treatment timepoints |
| `02-comparisons_with_healthy.ipynb` | MM vs. healthy donor protein levels |
| `03-pathway_analyses.ipynb` | Pathway enrichment on differential proteins |

## Notes

- Daratumumab-treated subjects are excluded from all objects.
- Healthy donors (BR1/BR2) are included in the Plasma object only.
- Only NDMM treatment and flu-series timepoints are retained.

## Environment

R

## Contributors

- Preliminary Olink analyses by [Love Tatting](mailto:love.tatting@alleninstitute.org)
- Final Olink analyses by Imran McGrath & Samir Rachid Zaim
