# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(data.table)
  library(stringr)
  library(ComplexHeatmap)
  library(circlize)
  library(glue)
  library(tidyverse)
})


#' Plot a gene-level heatmap for selected pathways across cell types
#'
#'
#' @param tissue tissue, either bmmc or pbmc
#' @param comp comparison to plot, ie PreTx_vs_Healthy
#' @param pathways pathways to include in heatmap
#' @param celltypes celltypes to include in heatmap
#' @param celltype_level celltype level, either "l2" or "l3"
#' @param title optional plot title; defaults to comp if not provided
#' @param low_col color for the low end of the heatmap scale (default: "blue")
#' @param mid_col color for the midpoint (zero) of the heatmap scale (default: "white")
#' @param high_col color for the high end of the heatmap scale (default: "red")
#'
#' @return Heatmap displaying log2FC values
#'   for the top leading-edge genes, with asterisks marking significant genes
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
    low_col = "blue",
    mid_col = "white",
    high_col = "red") {
  # use comp as title if not provided
  plot_title <- if (!is.null(title)) title else comp

  pathway_df <- fread(glue(
    "{input_dir}/{celltype_level}-fgsea_all_results/fgsea_all_results_{tissue}_{comp}.csv"
  )) %>%
    filter(pathway %in% pathways) %>%
    mutate(
      pathway = gsub("[^A-Za-z0-9_]", "_", pathway),
      pathway_db = paste(pathway, gmt_db, sep = "_")
    ) %>%
    filter(pathway_db != "OXIDATIVE_PHOSPHORYLATION_KEGG_2021_Human.gmt")
  pathways_clean <- pathway_df$pathway %>% unique()

  deseq <- fread(glue(
    "{input_dir}/{celltype_level}-deseq2_results/deseq2_results_{tissue}_{comp}.csv"
  )) %>%
    filter(!is.na(padj)) %>%
    mutate(sig = Qvalue < 0.05)

  celltypes <- intersect(celltypes, pathway_df$celltype)
  df_long <- pathway_df %>%
    filter(celltype %in% celltypes) %>%
    separate_rows(leadingEdge, sep = ",\\s*") %>%
    distinct(pathway, gene = leadingEdge)

  # get top leading edge genes per pathway
  final_genes <- lapply(pathways_clean, function(p) {
    pathway_genes <- df_long %>%
      filter(pathway == p) %>%
      pull(gene)

    deseq %>%
      filter(gene %in% pathway_genes & celltype %in% celltypes) %>%
      group_by(gene) %>%
      summarise(median_abs_lfc = median(abs(log2FoldChange))) %>%
      dplyr::slice_max(median_abs_lfc, n = 15) %>%
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


  # create deseq and pval matrices
  deseq_mat <- gene_pathway_matrix %>%
    left_join(deseq, by = "gene") %>%
    filter(celltype %in% celltypes & gene %in% final_genes) %>%
    select(gene, celltype, log2FoldChange) %>%
    pivot_wider(names_from = celltype, values_from = log2FoldChange) %>%
    tibble::column_to_rownames("gene") %>%
    as.matrix() %>%
    t()
  deseq_mat[is.na(deseq_mat)] <- 0
  deseq_mat <- deseq_mat[celltypes, ]

  pval_mat <- gene_pathway_matrix %>%
    left_join(deseq, by = "gene") %>%
    filter(celltype %in% celltypes & gene %in% final_genes) %>%
    select(gene, celltype, Qvalue) %>%
    pivot_wider(names_from = celltype, values_from = Qvalue) %>%
    tibble::column_to_rownames("gene") %>%
    as.matrix() %>%
    t()
  pval_mat <- pval_mat[celltypes, ]

  # create significance matrix for deseq values
  sig_mat <- ifelse(is.na(pval_mat), "",
    ifelse(pval_mat < 0.05, "*", "")
  )

  min_val <- min(deseq_mat, na.rm = TRUE)
  max_val <- max(deseq_mat, na.rm = TRUE)

  # fix white at 0; colors are user-configurable
  abs_max <- max(abs(c(min_val, max_val)), na.rm = TRUE)
  col_fun <- colorRamp2(
    c(-abs_max, 0, abs_max),
    c(low_col, mid_col, high_col)
  )

  cell_fun <- function(j, i, x, y, width, height, fill) {
    # use gray90 for NA values, otherwise use the color scale
    cell_fill <- if (is.na(pval_mat[i, j])) "gray90" else col_fun(deseq_mat[i, j])

    # base heatmap color with gray gridlines
    grid.rect(
      x = x, y = y, width = width, height = height,
      gp = gpar(col = "black", lwd = 0.5, fill = cell_fill)
    )

    # add dot if significant
    if (sig_mat[i, j] == "*") {
      grid.text("*", x = x, y = y, gp = gpar(fontsize = 15, col = "black"))
    }
  }
  rownames(gene_pathway_matrix) <- gene_pathway_matrix$gene

  # subset to genes actually present in deseq_mat
  gene_pathway_matrix_mat_sub <- gene_pathway_matrix[colnames(deseq_mat), , drop = FALSE] %>% select(!gene)

  # convert TRUE/FALSE to a factor for coloring
  gene_pathway_factors <- lapply(
    as.data.frame(gene_pathway_matrix_mat_sub),
    function(x) factor(x, levels = c(FALSE, TRUE))
  )

  uniform_colors <- list()
  for (p in colnames(gene_pathway_matrix_mat_sub)) {
    uniform_colors[[p]] <- c("FALSE" = "white", "TRUE" = "#77b8a1")
  }

  # cell_fun for annotation heatmaps to add gray gridlines
  anno_cell_fun <- function(j, i, x, y, width, height, fill) {
    grid.rect(
      x = x, y = y, width = width, height = height,
      gp = gpar(col = "black", lwd = 0.75, fill = fill)
    )
  }

  # create column annotation
  col_ha <- HeatmapAnnotation(
    df = as.data.frame(gene_pathway_factors),
    col = uniform_colors,
    annotation_name_side = "left",
    show_legend = FALSE,
    simple_anno_size_adjust = TRUE,
    gp = gpar(col = "black", lwd = 0.75)
  )

  Heatmap(
    deseq_mat,
    name = "Log2FC",
    cell_fun = cell_fun,
    column_title = plot_title,
    cluster_rows = F,
    cluster_columns = F,
    top_annotation = col_ha
  )
}
