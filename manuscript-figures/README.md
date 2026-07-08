# manuscript-figures/

Publication figure notebooks — one notebook per panel. Notebooks read pre-computed inputs from `data/figure-inputs/` and write to `data/figure-outputs/`.

## Structure

```
figure-1/ … figure-7/      # main figure panels
supplementary-1/ … 7/      # supplementary panels
utilities/                 # shared R plotting functions and environment
```

Batch-run with `run_py_nbs.sh` (Python panels) and `run_r_nbs.sh` (R panels).

## Environment

R (primary) — `utilities/environment/r-env.yml` (setup: `utilities/environment/instructions.txt`). Some panels use the scRNA Python environment.

## Shared plotting functions (`utilities/functions/`)

`label-and-color-maps.R` (cell-type labels, color palettes, orderings) is the base — source it at the top of any notebook:

```r
source('../utilities/functions/label-and-color-maps.R')
```

The remaining files provide frequency, DEG, pathway, and gene-expression plot layouts, plus Olink helpers (`helper_functions_olink.R`).
