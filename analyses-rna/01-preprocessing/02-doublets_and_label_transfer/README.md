## 02-doublets_and_label_transfer
#### Code by [Aishwarya Chander](mailto:aishwarya.chander@alleninstitute.org)

### Files

- `01-predict_doublets.py` — identifies doublets using `Scrublet`
- `02-predict_labels.py` — predicts cell-type labels using CellTypist against 3 references:
  - Healthy PBMC — AIFI Immune Health reference
  - Healthy BMMC — Flex-kit [BMMC reference](https://explore.allenimmunology.org/explore/e4fdcc96-ed2c-485a-a121-5ea91192355e)
  - External BMMC — [Bone Marrow Atlas](https://www.cell.com/cell/fulltext/S0092-8674(24)00408-2)

### Environment

`ndmm-scrna` (Python): `scrublet`, `celltypist`, `scanpy`, `pandas`. The `.py` files are percent-format notebooks; ensure sufficient memory for large single-cell datasets.
