# Multiple myeloma and therapy reshape the bone marrow niche to durably constrain immune reconstitution and vaccine responsiveness
## Data analyses and figures

- **Cohort:** NDMM patients sampled at Pre-Treatment, End of Induction, and Post-Transplant (60d, 90d, 1y, 2y), with matched healthy donors and flu-vaccination subcohorts.
- **Compartments:** BMMC (bone marrow) and PBMC (peripheral blood).

## Repository structure

```
analyses-rna/         # scRNA-seq: preprocessing, object building, subset & pseudobulk analyses
analyses-flow/        # Flow cytometry composition analysis
analyses-olink/       # Olink proteomics (plasma + bone marrow)
analyses-msd/         # MSD cytokine assay cleaning (vaccine responses)
manuscript-figures/   # Figure-panel notebooks
data/                 # All inputs and outputs, organized by modality (not tracked in git)
setup/                # Environment and data-staging helpers
```

## Data

Processed data and results are deposited on Zenodo; raw sequencing data are in GEO. Download into `data/` (see `setup/`); each `analyses-*` directory reads and writes under `data/<modality>/`.

## Primary comparisons

- **Longitudinal:** Pre-Treatment → End Induction → Post-Transplant (60d, 90d, 1y, 2y)
- **Cross-sectional:** MM vs. healthy donors
- **Responder groups:** responders vs. non-responders
- **Vaccine response:** flu and COVID series

## Team

| Name | Role | Modality |
|---|---|---|
| Aishwarya Chander | Analyst | scRNA-seq |
| Medbh Dillon | Research Associate | Flow |
| Palak Genge | Analyst | scRNA-seq |
| Nicolas Moss | Analyst | Figure support |
| Kevin Lee | Analyst | Flow |
| Melinda Angus-Hill | Senior Scientist | Oncology |
| Ziyuan He | Senior Scientist | Advisor |
| *Imran McGrath | Analyst | Olink |
| *Samir Rachid Zaim | Senior Scientist | Biostatistics |

*These contributors are Allen Institute Alumnis. 