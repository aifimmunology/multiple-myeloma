# analyses-flow/

Flow cytometry composition analysis of PBMC from the NDMM cohort. Quantifies T, B, and myeloid populations across treatment timepoints and the flu-vaccination series, compared to healthy donors and between responder groups.

## Scripts

| File | Description |
|---|---|
| `01_qc_processing_adjusted.ipynb` | Load and clean T/B/myeloid panels; filter viable cells, standardize timepoints, remove duplicates |
| `02-ndmm-vs-healthy-comparison.ipynb` | MM vs. healthy donors, cross-sectional across timepoints |
| `03-vrd-vs-healthy.ipynb` | VRd-treated vs. healthy |
| `04-dvrd-vs-healthy.ipynb` | DVRd-treated vs. healthy |

## Methods

- **Longitudinal:** CLR-transformed cell frequencies → paired tests across flu timepoints (SA → Day 0 → Day 7 → Day 90).
- **Cross-sectional:** CLR + linear regression controlling for sex and response — `CLR(frequency) ~ Sex + Response_Status`.
- **Timepoints compared:** Healthy vs Pre-Trx, EI, ASCT60d, ASCT90d, ASCT1y; Pre-Trx vs PI2C, EI; EI vs ASCT60d, ASCT90d, ASCT1y, ASCT2y.

## Environment

R (`r_env.sh`)

## Data

Inputs (sample key, cohort metadata, cleaned panels) and outputs (frequency tables, plots, stats) under `data/flow/`.

## Contributors

- Flow data annotated by [Medbh Dillon](mailto:medbh.dillon@alleninstitute.org)
- Preliminary analyses by [Kevin Lee](mailto:kevin.lee@alleninstitute.org)
- Final analyses by [Samir Rachid Zaim](mailto:samir.rachidzaim@alleninstitute.org)
