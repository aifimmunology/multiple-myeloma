# Author: Aishwarya Chander
# Organization: Allen Institute for Immunology

# Volcano plots for differential expression results

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggrepel)
    library(grid)
})

source("../utilities/functions/label-and-color-maps.R")

#' Create faceted volcano plots for DEG results
#' @param input_deg_res Path to DEG results CSV
#' @param selected_celltypes Cell types to include
#' @param timepoint_1 First timepoint (positive fold change direction)
#' @param timepoint_2 Second timepoint (negative fold change direction)
#' @param ncols Number of columns in facet grid
#' @param file_save_name Optional path to save plot
#' @param plot_width Plot width (auto-calculated if NULL)
#' @param plot_height Plot height (auto-calculated if NULL)
#' @param max_overlaps Maximum number of overlaps allowed for ggrepel labels
#' @param save_dpi DPI for saved plot
#' @return ggplot object
plot_volcano_publication <- function(input_deg_res,
                                     selected_celltypes,
                                     timepoint_1 = "PreTx",
                                     timepoint_2 = "Healthy",
                                     ncols = 5,
                                     file_save_name = NULL,
                                     title = NULL,
                                     plot_width = NULL,
                                     plot_height = NULL,
                                     max_overlaps = 20,
                                     save_dpi = 200) {

    if (!file.exists(input_deg_res)) stop("File not found: ", input_deg_res)
    res <- fread(input_deg_res)

    required_cols <- c("celltype", "Qvalue", "log2FoldChange", "pvalue", "gene")
    missing_cols <- setdiff(required_cols, colnames(res))
    if (length(missing_cols)) {
        stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
    }

    color_tp_1 <- get_timepoint_colors(timepoint_1)
    color_tp_2 <- get_timepoint_colors(timepoint_2)

    if (is.na(color_tp_1)) stop("Color not found for timepoint_1: ", timepoint_1)
    if (is.na(color_tp_2)) stop("Color not found for timepoint_2: ", timepoint_2)

    res <- res[!is.na(Qvalue) & !is.na(log2FoldChange) & celltype %in% selected_celltypes]
    if (nrow(res) == 0) stop("No rows after filtering for selected cell types.")

    res[, celltype := factor(as.character(celltype), levels = selected_celltypes)]
    res[, significance := "Not Sig"]
    res[Qvalue < 0.1 & log2FoldChange > 0.5,  significance := paste("Up in", timepoint_1)]
    res[Qvalue < 0.1 & log2FoldChange < -0.5, significance := paste("Up in", timepoint_2)]

    sig_levels <- c(paste("Up in", timepoint_2), "Not Sig", paste("Up in", timepoint_1))
    res[, significance := factor(significance, levels = sig_levels)]
    manual_colors <- setNames(c(color_tp_2, "gray70", color_tp_1), sig_levels)

    n_celltypes <- length(selected_celltypes)
    nrows <- ceiling(n_celltypes / ncols)
    if (is.null(plot_width)) plot_width <- 3 * ncols
    if (is.null(plot_height)) plot_height <- 4 * nrows

    p <- ggplot(res, aes(x = log2FoldChange, y = -log10(pvalue),
                          color = significance, label = gene)) +
        geom_point(size = 1, alpha = 0.8) +
        geom_text_repel(
            data = res[Qvalue < 0.1 & abs(log2FoldChange) > 0.5],
            force = 1.5, max.overlaps = max_overlaps, size = 3.5, show.legend = FALSE
        ) +
        geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed",
                   color = "grey50", linewidth = 0.5) +
        facet_wrap(~ celltype, ncol = ncols, scales = "free") +
        scale_color_manual(values = manual_colors) +
        labs(
            title = title,
            x = expression(Log[2]~Fold~Change),
            y = expression(-Log[10]~"(p-value)"),
            color = NULL
        ) +
        theme_classic() +
        theme(
                plot.title = element_text(size = 17, face = "bold"),
            strip.text = element_text(size = 15, hjust = 0, face = "bold"),
            strip.background = element_blank(),
            axis.title = element_text(size = 15, face = "bold"),
            axis.text = element_text(size = 12),
            legend.position = "bottom",
            legend.text = element_text(size = 11),
            panel.spacing = unit(0.5, "cm")
        ) +
        guides(color = guide_legend(override.aes = list(size = 4)))

    if (!is.null(file_save_name)) {
        ggsave(filename = file_save_name, plot = p,
               width = plot_width, height = plot_height, dpi = save_dpi)
    }

    return(p)
}
