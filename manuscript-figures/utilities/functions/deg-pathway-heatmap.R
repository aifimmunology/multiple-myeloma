# Author: Nicholas Moss
# Edited by: Aishwarya Chander
# Organization: Allen Institute for Immunology

# Combined DEG heatmap + pathway annotation figure
# Refactored from rxiv-figs/r-deg_plots.ipynb

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

source("../utilities/functions/label-and-color-maps.R")

#' Plot a combined DEG heatmap with fgsea pathway annotation
#'
#' Top panel: binary pathway membership grid (leading-edge genes x pathways).
#' Bottom panel: log2FC heatmap (genes x cell types), grey = untested,
#'   white dots = padj < sig_padj_cut.
#' Gene order uses a pathway-aware staircase strategy; pathways are
#' hierarchically clustered on binary membership.
#'
#' @param deg_file        Path to DESeq2 results CSV (requires: gene, celltype,
#'                        log2FoldChange, and optionally padj)
#' @param fgsea_file      Path to fgsea results CSV (requires: pathway, celltype,
#'                        and a leading-edge column matching "leading.?edge")
#' @param celltypes       Character vector of cell types to include
#' @param lfc_cut         Minimum |log2FoldChange| for a gene to pass the filter
#'                        (applied per gene across any of the requested celltypes)
#' @param padj_cut        If not NULL, genes must also have padj < padj_cut in at
#'                        least one celltype to pass the filter (default NULL)
#' @param cap             LFC values are capped to [-cap, cap] for the color scale
#' @param top_n_pw        Number of top pathways to show (ranked by leading-edge
#'                        overlap with the filtered gene set); ignored when
#'                        \code{pathways} is supplied
#' @param pathways        Optional character vector of pathway names to display.
#'                        When provided, bypasses \code{top_n_pw} auto-selection
#'                        and uses these pathways directly. Order is preserved.
#' @param sig_padj_cut    padj threshold for the white significance dots
#'                        (default 0.05; ignored if padj column absent)
#' @param timepoint_1     Numerator timepoint name (positive LFC = up here);
#'                        looked up in tp_color_map for color_high
#' @param timepoint_2     Denominator timepoint name (negative LFC = up here);
#'                        looked up in tp_color_map for color_low
#' @param fill_label      Legend title for the LFC color bar; auto-generated as
#'                        "log2FC\\n(t1 vs t2)" if NULL
#' @param untested_fill   Fill color for untested gene x celltype tiles
#' @param plot_width      repr plot width (passed to options())
#' @param plot_height     repr plot height (passed to options())
#' @param panel_heights   Relative heights of annotation and heatmap panels
#' @param base_size       Base font size for both panels
#' @param x_text_size     Font size for x-axis gene labels
#' @param y_text_size     Font size for y-axis celltype / pathway labels
#' @return patchwork ggplot object

plot_deg_pathway_heatmap <- function(
    deg_file,
    fgsea_file,
    celltypes,
    lfc_cut = 1,
    padj_cut = NULL,
    cap = 2,
    top_n_pw = 10,
    pathways = NULL,
    sig_padj_cut = 0.05,
    timepoint_1 = NULL,
    timepoint_2 = NULL,
    fill_label = NULL,
    untested_fill = "grey85",
    plot_width = 20,
    plot_height = 6,
    panel_heights = c(1, 1),
    base_size = 10,
    x_text_size = 12,
    y_text_size = 12) {
  options(repr.plot.width = plot_width, repr.plot.height = plot_height)

  # 0. Resolve colors from tp_color_map
  color_high <- if (!is.null(timepoint_1) && timepoint_1 %in% names(tp_color_map)) {
    tp_color_map[[timepoint_1]]
  } else {
    "red"
  }

  color_low <- if (!is.null(timepoint_2) && timepoint_2 %in% names(tp_color_map)) {
    tp_color_map[[timepoint_2]]
  } else {
    "blue"
  }

  if (is.null(fill_label)) {
    fill_label <- if (!is.null(timepoint_1) && !is.null(timepoint_2)) {
      paste0("log2FC\n(", timepoint_1, " vs ", timepoint_2, ")")
    } else {
      "log2FC"
    }
  }

  # 1. Load DEG data
  deg <- read_csv(deg_file, show_col_types = FALSE)
  has_padj <- "padj" %in% colnames(deg)

  deg_sub <- deg %>% filter(celltype %in% celltypes)

  # 2. Filter genes
  passing_genes <- deg_sub %>%
    group_by(gene) %>%
    summarize(
      n = sum(
        abs(log2FoldChange) >= lfc_cut &
          if (!is.null(padj_cut) && has_padj) (!is.na(padj) & padj < padj_cut) else TRUE,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    filter(n >= 1) %>%
    pull(gene)

  if (length(passing_genes) == 0) stop("No genes passed the filter — try lowering lfc_cut or padj_cut.")

  deg_sub <- deg_sub %>% filter(gene %in% passing_genes)

  # Full gene x celltype grid — NA = untested, not imputed
  deg_full <- expand_grid(gene = unique(passing_genes), celltype = celltypes) %>%
    left_join(
      deg_sub %>% select(gene, celltype, log2FoldChange, any_of("padj")),
      by = c("gene", "celltype")
    ) %>%
    mutate(
      tested = !is.na(log2FoldChange),
      padj   = if (has_padj) padj else NA_real_
    )

  # 3. Parse fgsea leading edges
  parse_le <- function(x) {
    if (is.na(x)) {
      return(character(0))
    }
    s <- gsub("c\\(|\\)|\\[|\\]|\"|'", "", as.character(x))
    strsplit(trimws(s), "[,;\\s]+", perl = TRUE)[[1]]
  }

  fg <- read_csv(fgsea_file, show_col_types = FALSE) %>%
    filter(celltype %in% celltypes)

  le_col <- grep("leading.?edge", names(fg), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(le_col)) stop("No leading-edge column found in fgsea results.")

  pw_genes <- fg %>%
    transmute(pathway, le = .data[[le_col]]) %>%
    mutate(gene = lapply(le, parse_le)) %>%
    select(pathway, gene) %>%
    unnest(gene) %>%
    filter(gene %in% passing_genes) %>%
    distinct()

  # 4. Select top pathways
  if (!is.null(pathways)) {
    missing_pw <- setdiff(pathways, unique(fg$pathway))
    if (length(missing_pw)) {
      warning(
        "Pathways not found in fgsea results and will be skipped: ",
        paste(missing_pw, collapse = ", ")
      )
    }
    top_pws <- intersect(pathways, unique(pw_genes$pathway))
  } else {
    top_pws <- pw_genes %>%
      count(pathway, sort = TRUE) %>%
      slice_head(n = top_n_pw) %>%
      pull(pathway)
  }

  if (length(top_pws) == 0) stop("No pathways found overlapping the filtered gene set.")

  # 5. Cluster pathways
  pw_df <- pw_genes %>%
    filter(pathway %in% top_pws) %>%
    mutate(val = 1L) %>%
    pivot_wider(names_from = gene, values_from = val, values_fill = 0L) %>%
    tibble::column_to_rownames("pathway")

  missing_cols <- setdiff(passing_genes, colnames(pw_df))
  if (length(missing_cols)) pw_df[, missing_cols] <- 0L

  pw_mat <- as.matrix(pw_df)
  pw_hc <- hclust(
    dist(pw_mat[, intersect(passing_genes, colnames(pw_mat))], method = "binary"),
    method = "ward.D2"
  )
  pw_order <- rownames(pw_mat)[pw_hc$order]

  pathway_symbols <- c(  "●", "▲", "▼", "■", "◆", "★", "✚", "✦", "✖", "✱", "⬟", "⬢", "♦", "♣", "♥")
  pathway_symbol_map <- setNames(
    sample(
      pathway_symbols,
      length(pw_order),
      replace = length(pw_order) > length(pathway_symbols)
    ),
    pw_order
  )
  pathway_labels <- setNames(
    paste(pathway_symbol_map[pw_order], pw_order),
    pw_order
  )

  # 6. Order genes (pathway-aware staircase)
  pw_gene_membership <- pw_genes %>%
    filter(pathway %in% top_pws) %>%
    mutate(pw_rank = match(pathway, pw_order)) %>%
    group_by(gene) %>%
    summarize(top_pw = min(pw_rank), .groups = "drop")

  all_genes_ranked <- tibble(gene = passing_genes) %>%
    left_join(pw_gene_membership, by = "gene") %>%
    mutate(top_pw = replace_na(top_pw, length(pw_order) + 1L))

  # Use 0 for NAs only for clustering, not for display
  lfc_mat <- deg_full %>%
    mutate(lfc_for_clust = replace_na(log2FoldChange, 0)) %>%
    select(gene, celltype, lfc_for_clust) %>%
    pivot_wider(names_from = gene, values_from = lfc_for_clust) %>%
    tibble::column_to_rownames("celltype") %>%
    as.matrix()
  lfc_mat <- lfc_mat[celltypes, , drop = FALSE]

  gene_order <- all_genes_ranked %>%
    group_by(top_pw) %>%
    group_modify(~ {
      genes <- .x$gene
      if (length(genes) <= 2) {
        return(.x %>% mutate(sub_order = seq_along(gene)))
      }
      sub_mat <- t(lfc_mat[, genes, drop = FALSE])
      hc <- hclust(dist(sub_mat), method = "ward.D2")
      .x %>% mutate(sub_order = match(gene, genes[hc$order]))
    }) %>%
    ungroup() %>%
    arrange(top_pw, sub_order) %>%
    pull(gene)

  # 7. Build plot data
  plot_df <- deg_full %>%
    mutate(
      gene     = factor(gene, levels = gene_order),
      celltype = factor(celltype, levels = celltypes),
      lfc_cap  = pmax(pmin(log2FoldChange, cap), -cap),
      sig      = has_padj & !is.na(padj) & padj < sig_padj_cut
    )

  ann_df <- pw_genes %>%
    filter(pathway %in% top_pws) %>%
    mutate(
      gene    = factor(gene, levels = gene_order),
      pathway = factor(pathway, levels = pw_order)
    )

  # 8. Pathway annotation panel (top)
  ann_grid <- expand_grid(
    gene    = factor(gene_order, levels = gene_order),
    pathway = factor(pw_order, levels = pw_order)
  )

  p_ann <- ggplot() +
    geom_tile(
      data = ann_grid, aes(x = gene, y = pathway),
      fill = "white", color = "black", linewidth = 0.1
    ) +
    geom_tile(
      data = ann_df, aes(x = gene, y = pathway),
      fill = "#edb779", color = "black", linewidth = 0.4
    ) +
    scale_y_discrete(limits = pw_order, labels = pathway_labels, drop = FALSE) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = base_size) +
    theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line    = element_blank(),
      axis.text.y  = element_text(size = y_text_size),
      plot.margin  = margin(5, 5.5, 0, 5.5)
    )

  # 9. LFC heatmap panel (bottom)
  not_tested_df <- plot_df %>% filter(!tested)
  tested_df <- plot_df %>% filter(tested)

  p_heat <- ggplot() +
    geom_tile(
      data = not_tested_df,
      aes(x = gene, y = celltype),
      fill = untested_fill, color = "black", linewidth = 0.1
    ) +
    geom_tile(
      data = tested_df,
      aes(x = gene, y = celltype, fill = lfc_cap),
      color = "black", linewidth = 0.4
    ) +
    {
      if (has_padj) {
        geom_point(
          data = filter(tested_df, sig),
          aes(x = gene, y = celltype),
          size = 2, color = "black", shape = 8, stroke = 0.5
        )
      }
    } +
    scale_fill_gradient2(
      low      = color_low,
      mid      = "white",
      high     = color_high,
      midpoint = 0,
      limits   = c(-cap, cap),
      na.value = untested_fill,
      name     = fill_label
    ) +
    labs(x = "DEGs (leading-edge filtered)", y = NULL, size = 5) +  
    theme_minimal(base_size = base_size) +
    theme(
      axis.line = element_blank(),
      panel.grid = element_blank(),
      axis.text.y = element_text(size = y_text_size, color = "black"),
      axis.text.x = element_text(
        size = x_text_size, angle = 90,
        vjust = 0.5, hjust = 1, color = "black"
      ),
      axis.title.x = element_text(size = 14),
      plot.margin = margin(0, 5.5, 5, 5.5)
    )

  # 10. Combine
  (p_ann / p_heat) + plot_layout(heights = panel_heights)
}
