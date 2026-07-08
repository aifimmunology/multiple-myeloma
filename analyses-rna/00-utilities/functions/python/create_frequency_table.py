# This code processes single-cell metadata to compute cell type counts and frequencies per sample, applies CLR (centered log-ratio) transformation for compositional data analysis, and merges sample-level metadata. Designed for analyzing cell type proportions across samples/conditions.

import pandas as pd
import numpy as np
import re

def clr_transform(x):
    x = np.array(x)
    if len(x) == 0 or np.any(x <= 0):
        return np.full_like(x, np.nan, dtype=np.float64)
    gm = np.exp(np.mean(np.log(x)))
    return np.log(x / gm)

def clean_dtypes(obj_meta: pd.DataFrame) -> pd.DataFrame:
    # Convert ints/categoricals to object/float to avoid groupby dtype surprises
    cast_map = {}
    for col in obj_meta.columns:
        if pd.api.types.is_integer_dtype(obj_meta[col]):
            cast_map[col] = "float"
        elif pd.api.types.is_categorical_dtype(obj_meta[col]):
            cast_map[col] = "object"
    return obj_meta.astype(cast_map) if cast_map else obj_meta

def compute_counts_with_clr(obj_meta: pd.DataFrame, groupby_col: str, celltype_col: str) -> pd.DataFrame:
    # Capture l1, l1.5, l2, l3...
    m = re.search(r'(l\d(?:\.\d)?)', celltype_col)
    if not m:
        raise ValueError(f"Could not extract level id from {celltype_col}")
    level_id = m.group(1)

    # Ensure full combinations are preserved
    obj_meta[celltype_col] = pd.Categorical(obj_meta[celltype_col])
    obj_meta[groupby_col] = pd.Categorical(obj_meta[groupby_col])

    counts = (
        obj_meta.groupby([celltype_col, groupby_col], observed=False)
        .size().reset_index(name="counts")
    )

    # Fill in missing sample×cell-type with 0
    full = pd.MultiIndex.from_product(
        [obj_meta[celltype_col].cat.categories, obj_meta[groupby_col].cat.categories],
        names=[celltype_col, groupby_col]
    ).to_frame(index=False)

    counts = (full.merge(counts, on=[celltype_col, groupby_col], how="left")
                   .fillna({"counts": 0}))

    # Raw frequency
    totals = counts.groupby(groupby_col, observed=False)["counts"].sum().reset_index(name="total_counts")
    counts = counts.merge(totals, on=groupby_col, how="left")
    counts["raw_frequency"] = counts["counts"] / counts["total_counts"].replace({0: np.nan})

    # Pseudocount frequency --> Pseudocount of 1 prevents log(0) in CLR transformation for zero-count cell types
    counts["pseudo_counts"] = counts["counts"] + 1
    pseudo_totals = (counts.groupby(groupby_col, observed=False)["pseudo_counts"]
                           .sum().reset_index(name="pseudo_total_counts"))
    counts = counts.merge(pseudo_totals, on=groupby_col, how="left")
    counts["pseudo_frequency"] = counts["pseudo_counts"] / counts["pseudo_total_counts"]

    # CLR across cell types within each sample
    counts["cell_type_clr"] = counts.groupby(groupby_col, observed=False)["pseudo_frequency"]\
                                    .transform(lambda x: clr_transform(x.values))

    counts["level_id"] = level_id
    return counts

def extract_and_merge_metadata(
    obj_meta: pd.DataFrame,
    counts_df: pd.DataFrame,
    desired_cols: list,
    groupby_col: str,
    visit_col: str | None = None,
    visit_order: list | None = None,
) -> pd.DataFrame:
    # Ensure groupby_col is available to merge back
    if groupby_col not in desired_cols:
        desired_cols = [groupby_col] + [c for c in desired_cols if c != groupby_col]

    existing = [c for c in desired_cols if c in obj_meta.columns]
    md = obj_meta[existing].drop_duplicates().reset_index(drop=True).astype(object)

    # Add any missing desired cols as "NA"
    for c in desired_cols:
        if c not in md.columns:
            md[c] = "NA"

    if visit_col and visit_col in md.columns and visit_order:
        md[visit_col] = pd.Categorical(md[visit_col], categories=visit_order, ordered=True)

    merged = counts_df.merge(md, on=groupby_col, how="left")
    return merged

def process_single_cell_metadata(
    adata,
    groupby_col: str,
    celltype_col: str,
    desired_cols: list,
    visit_col: str | None = None,
    visit_order: list | None = None,
) -> pd.DataFrame:
    obj_meta = clean_dtypes(adata.obs.copy())
    counts_df = compute_counts_with_clr(obj_meta, groupby_col, celltype_col)
    merged = extract_and_merge_metadata(
        obj_meta=obj_meta,
        counts_df=counts_df,
        desired_cols=desired_cols,
        groupby_col=groupby_col,
        visit_col=visit_col,
        visit_order=visit_order,
    )
    return merged