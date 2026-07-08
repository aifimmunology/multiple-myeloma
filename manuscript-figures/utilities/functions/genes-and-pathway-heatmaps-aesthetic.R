# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(stringr)
  library(glue)
  library(tidyverse)
})

source("../utilities/functions/label-and-color-maps.R")

#' Plot a gene-level heatmap for selected pathways across cell types
#'
#'
#' @param tissue tissue, either bmmc or pbmc
#' @param comp comparison to plot, ie PreTx_vs_Healthy
#' @param pathways pathways to include in heatmap
#' @param celltypes celltypes to include in heatmap
#' @param celltype_level celltype level, either "l2" or "l3"
#' @param title optional plot title; defaults to comp if not provided
#' @param timepoint_1     Numerator timepoint name (positive LFC = up here);
#'                        looked up in tp_color_map for color_high
#' @param timepoint_2     Denominator timepoint name (negative LFC = up here);
#'                        looked up in tp_color_map for color_low
#' @param fill_label      Legend title for the LFC color bar; auto-generated as
#'                        "log2FC\n(t1 vs t2)" if NULL
#' @param cap LFC values are capped to [-cap, cap] for the color scale
#' @param sig_padj_cut padj threshold for significance dots (default 0.05)
#' @param untested_fill Fill color for untested gene x celltype tiles
#' @param plot_width repr plot width (passed to options())
#' @param plot_height repr plot height (passed to options())
#' @param panel_heights Relative heights of annotation and heatmap panels
#' @param base_size Base font size for both panels
#' @param x_text_size Font size for x-axis gene labels
#' @param y_text_size Font size for y-axis celltype / pathway labels
#' @param top_n_genes Number of top leading-edge genes to keep per pathway
#'
#' @return patchwork ggplot object displaying log2FC values
#'   for the top leading-edge genes, with white dots marking significant genes
#'   and a pathway-membership annotation along the top.
#'

plot_pathway_genes_heatmap <- function(
    tissue = NULL,
    comp = NULL,
    pathways = NULL,
    celltypes = NULL,
    celltype_level = "l3",
    input_dir = "02-inputs",
    title = NULL,
    timepoint_1 = NULL,
    timepoint_2 = NULL,
    fill_label = NULL,
    cap = 2,
    sig_padj_cut = 0.05,
    sig_color = "black", 
    untested_fill = "grey85",
    plot_width = 20,
    plot_height = 6,
    panel_heights = c(1, 1),
    base_size = 10,
    x_text_size = 12,
    y_text_size = 12,
    top_n_genes = 15) {
  options(repr.plot.width = plot_width, repr.plot.height = plot_height)

  # ── 0. Resolve colors from tp_color_map ───────────────────────────────────
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

  # use comp as title if not provided
  plot_title <- if (!is.null(title)) title else comp
  pathways_clean_input <- gsub("[^A-Za-z0-9_-]", " ", pathways)

  pathway_df <- fread(glue(
    "{input_dir}/{celltype_level}-fgsea_all_results/fgsea_all_results_{tissue}_{comp}.csv"
  )) %>%
    mutate(
      pathway = gsub("[^A-Za-z0-9_-]", " ", pathway),
      pathway_db = paste(pathway, gmt_db, sep = "_")
    ) %>%
    filter(pathway %in% pathways_clean_input) %>%
    filter(pathway_db != "OXIDATIVE PHOSPHORYLATION_KEGG_2021_Human.gmt")
  pathways_clean <- pathway_df$pathway %>% unique()

  deseq <- fread(glue(
    "{input_dir}/{celltype_level}-deseq2_results/deseq2_results_{tissue}_{comp}.csv"
  )) %>%
    filter(!is.na(Qvalue)) %>%
    mutate(sig = Qvalue < 0.05)

  celltypes <- intersect(celltypes, pathway_df$celltype)
  leading_edge_long <- pathway_df %>%
    filter(celltype %in% celltypes, !is.na(leadingEdge)) %>%
    mutate(
      leadingEdge = as.character(leadingEdge),
      leadingEdge = str_remove(leadingEdge, "^c\\("),
      leadingEdge = str_remove(leadingEdge, "\\)$")
    ) %>%
    separate_rows(leadingEdge, sep = "\\s*,\\s*") %>%
    mutate(
      gene = leadingEdge %>%
        str_trim() %>%
        str_remove_all("^[\\\"']|[\\\"']$") %>%
        str_trim()
    ) %>%
    filter(gene != "") %>%
    distinct(pathway, celltype, gene)

  df_long <- leading_edge_long %>%
    distinct(pathway, gene)

  # get top leading edge genes per pathway
  final_genes <- lapply(pathways_clean, function(p) {
    pathway_genes <- df_long %>%
      filter(pathway == p) %>%
      pull(gene)

    deseq %>%
      filter(gene %in% pathway_genes & celltype %in% celltypes) %>%
      group_by(gene) %>%
      summarise(median_abs_lfc = median(abs(log2FoldChange))) %>%
        dplyr::slice_max(median_abs_lfc, n = top_n_genes) %>%
      pull(gene)
  }) %>%
    unlist() %>%
    unique()

  # filter for genes and create wide df
  gene_pathway_matrix <- df_long %>%
    filter(gene %in% final_genes) %>%
    mutate(present = TRUE) %>%
    pivot_wider(
      names_from = pathway,
      values_from = present,
      values_fill = FALSE
    ) %>%
    as.data.frame()

  # ── Build pathway annotation data ──────────────────────────────────────────
  pw_genes_long <- df_long %>%
    filter(gene %in% final_genes)

  # ── Cluster pathways (hierarchical, binary distance) ─────────────────────
  pw_df <- pw_genes_long %>%
    mutate(val = 1L) %>%
    pivot_wider(names_from = gene, values_from = val, values_fill = 0L) %>%
    tibble::column_to_rownames("pathway")

  missing_cols <- setdiff(final_genes, colnames(pw_df))
  if (length(missing_cols)) pw_df[, missing_cols] <- 0L

  pw_mat <- as.matrix(pw_df)
  pw_hc <- hclust(
    dist(pw_mat[, intersect(final_genes, colnames(pw_mat))], method = "binary"),
    method = "ward.D2"
  )
  pw_order <- rownames(pw_mat)[pw_hc$order]

  pathway_symbols <- c(  "●", "▲", "▼", "■", "◆", "★", "✚", "✦", "✖", "✱", "⬟", "♣", "♥")
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

  # ── Build full gene x celltype grid ────────────────────────────────────────
  deg_full <- expand_grid(gene = final_genes, celltype = celltypes) %>%
    left_join(
      deseq %>% select(gene, celltype, log2FoldChange, Qvalue),
      by = c("gene", "celltype")
    ) %>%
    mutate(
      tested  = !is.na(log2FoldChange),
      lfc_cap = pmax(pmin(log2FoldChange, cap), -cap),
      sig     = !is.na(Qvalue) & Qvalue < sig_padj_cut
    )

  # ── Order genes (pathway-aware staircase) ────────────────────────────────
  pw_gene_membership <- pw_genes_long %>%
    mutate(pw_rank = match(pathway, pw_order)) %>%
    group_by(gene) %>%
    summarize(top_pw = min(pw_rank), .groups = "drop")

  all_genes_ranked <- tibble(gene = final_genes) %>%
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

  # ── Build plot data ────────────────────────────────────────────────────────
  plot_df <- deg_full %>%
    mutate(
      gene     = factor(gene, levels = gene_order),
      celltype = factor(celltype, levels = celltypes)
    )

  ann_df <- leading_edge_long %>%
    filter(gene %in% gene_order, pathway %in% pw_order) %>%
    distinct(pathway, gene) %>%
    mutate(
      gene    = factor(gene, levels = gene_order),
      pathway = factor(pathway, levels = pw_order)
    )

  # ── Pathway annotation panel (top) ─────────────────────────────────────────
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
      fill = "#77b8a1", color = "black", linewidth = 0.4
    ) +
    scale_x_discrete(limits = gene_order, drop = FALSE) +
    scale_y_discrete(limits = pw_order, labels = pathway_labels, drop = FALSE) +
    labs(x = NULL, y = NULL, title = plot_title) +
    theme_classic(base_size = base_size) +
    theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line    = element_blank(),
      axis.text.y  = element_text(size = y_text_size, color = "black"),
      plot.margin  = margin(5, 5.5, 0, 5.5)
    )

  # ── LFC heatmap panel (bottom) ────────────────────────────────────────────
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
    geom_point(
      data = filter(tested_df, sig),
      aes(x = gene, y = celltype),
      size = 2, color = sig_color, shape = 8, stroke = 0.5
    ) +
    scale_fill_gradient2(
      low      = color_low,
      mid      = "white",
      high     = color_high,
      midpoint = 0,
      limits   = c(-cap, cap),
      na.value = untested_fill,
      name     = fill_label
    ) +
    scale_x_discrete(limits = gene_order, drop = FALSE) +
    scale_y_discrete(position = "right", drop = FALSE) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = base_size) +
    theme(
      axis.line = element_blank(),
      panel.grid = element_blank(),
      axis.text.y = element_text(size = y_text_size, color = "black", hjust = 0),
      axis.text.x = element_text(
        size = x_text_size, angle = 90,
        vjust = 0.5, hjust = 1, color = "black"
      ),
      axis.title.x = element_text(size = 14),
      plot.margin = margin(0, 5.5, 5, 5.5)
    )

  # ── Combine ────────────────────────────────────────────────────────────────
  p_ann / p_heat + plot_layout(heights = panel_heights)
}
