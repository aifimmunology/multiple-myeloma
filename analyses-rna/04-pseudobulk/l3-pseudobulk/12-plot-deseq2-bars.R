# ---
# jupyter:
#   jupytext:
#     formats: R:percent
#     text_representation:
#       extension: .R
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: R (seurat_v5)
#     language: R
#     name: seurat_v5
# ---

# %%
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

# %%
#' Count Differentially Expressed Genes (DEGs) by Term
#'
#' Counts the number of upregulated and downregulated DEGs for each term in a
#' differential expression results dataframe. Optionally returns a bar plot or
#' the summarized dataframe.
#'
#' @param deg_res Dataframe containing differential expression results.
#' @param p_col Character, the column name for p-values (default: `"qval"`).
#' @param p_cutoff Numeric, the p-value cutoff for significance (default:
#' `0.05`).
#' @param es_col Character, the column name for effect sizes (default:
#' `"estimate"`).
#' @param es_cutoff Numeric, the effect size cutoff for significance (default:
#' `0`).
#' @param term_col Character, the column name for terms (default: `"term"`).
#' @param title Character, the title of the plot (default: `""`).
#' @param return_df Logical, whether to return the summarized dataframe instead
#' of the plot (default: `FALSE`).
#' @param label_nudge Numeric, the amount to nudge text labels on the plot
#' (default: `50`).
#' @param base_size Numeric, the base font size for the plot (default: `16`).
#' @param label_size Numeric, the size of the text labels (default: `5`).
#' @param facet_col Character, the column name for faceting the plot (default:
#' `NULL`).
#' @param ... Additional arguments passed to `facet_wrap` if `facet_col` is
#' specified.
#'
#' @return A ggplot object showing the bar plot of DEG counts or a dataframe
#' with the summarized counts if `return_df = TRUE`.
#'
#' @examples
#' # Example usage:
#' deg_res <- data.frame(
#'   term = c("Term1", "Term1", "Term2", "Term2"),
#'   estimate = c(0.5, -0.3, 0.6, -0.2),
#'   qval = c(0.01, 0.2, 0.03, 0.1)
#' )
#' deg_counts(deg_res, p_col = "qval", es_col = "estimate", term_col = "term")
deg_counts <- function(
  deg_res,
  p_col = "qval",
  p_cutoff = 0.05,
  es_col = "estimate",
  es_cutoff = 0,
  term_col = "term",
  title = "",
  return_df = FALSE,
  label_nudge = 50,
  base_size = 16,
  label_size = 5,
  facet_col = NULL,
  ...
) {
  # handle globals
  direction <- n_direction <- NULL

  term_vals <- unique(deg_res[[term_col]])
  term_vals <- term_vals[!grepl("Intercept|__", term_vals)]
  df <- deg_res %>%
    dplyr::filter(get(p_col) <= p_cutoff, abs(get(es_col)) >= es_cutoff) %>%
    dplyr::filter(!grepl("Intercept|__", get(term_col))) %>%
    dplyr::mutate(direction = sign(get(es_col))) %>%
    dplyr::group_by(across(all_of(c(term_col, facet_col))), direction) %>%
    dplyr::tally() %>%
    dplyr::mutate(n_direction = direction * n)
  df[[term_col]] <- factor(df[[term_col]], levels = term_vals)
  if (return_df) {
    df
  } else {
    sub_title <- sprintf("%s \u2264 %s; abs(%s) \u2265 %s", p_col,
                         p_cutoff, es_col, es_cutoff)
    g <- ggplot(df, aes(!!as.name(term_col), n_direction,
                        fill = factor(direction))) +
      geom_col() +
      theme_bw(base_size = base_size) +
      ylab("N Differential Genes") +
      scale_fill_manual(
        values = c("blue4", "red2"), breaks = c(-1, 1),
        labels = c("Down", "Up"), name = "Direction"
      ) +
      geom_text(
        data = df %>% filter(direction == -1),
        aes(label = abs(n_direction)), vjust = 0.5, size = label_size,
        nudge_y = -1 * label_nudge
      ) +
      geom_text(
        data = df %>% filter(direction == 1), size = label_size,
        aes(label = abs(n_direction)), vjust = 0.5,
        nudge_y = label_nudge
      ) +
      ggtitle(title, subtitle = sub_title) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    if (!is.null(facet_col)) {
      ff <- as.formula(paste0("~", facet_col))
      g <- g + facet_wrap(ff, ...)
    }
    g
  }
}


# %%
# assumes your deg_counts() is already defined
batch_deg_counts <- function(
    files,
    out_dir = "deg_counts_plots",
    p_col = "padj",
    p_cutoff = 0.05,
    es_col = "log2FoldChange",
    es_cutoff = 0,
    term_col = "contrast",
    facet_col = NULL, # e.g., "celltype" or "contrast"
    width = 15, height = 6, dpi = 300,
    also_pdf = FALSE) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  for (f in files) {
    message("Processing: ", f)
    dt <- fread(f)

    # Basic column checks with clear errors
    need <- c(p_col, es_col, term_col)
    missing <- setdiff(need, names(dt))
    if (length(missing)) {
      warning(sprintf("Skipping %s; missing columns: %s", f, paste(missing, collapse = ", ")))
      next
    }
    if (!is.null(facet_col) && !facet_col %in% names(dt)) {
      warning(sprintf("Facet column '%s' not found in %s; proceeding without faceting.", facet_col, f))
      facet_col <- NULL
    }

    title <- str_replace(basename(f), "_deg\\.csv$", "") |> str_replace("\\.csv$", "")

    dt[[paste0(term_col, "_label")]] <- factor(
      stringr::str_to_title(gsub("_", " ", dt[[term_col]])),
      levels = stringr::str_to_title(gsub("_", " ", levels(factor(dt[[term_col]]))))
    )
    term_col_label <- paste0(term_col, "_label")

    g <- deg_counts(
      deg_res   = dt,
      p_col     = p_col,
      p_cutoff  = p_cutoff,
      es_col    = es_col,
      es_cutoff = es_cutoff,
      term_col  = term_col_label,
      title     = title,
      return_df = FALSE,
      facet_col = facet_col
    )

    png_file <- file.path(out_dir, paste0(title, "_deg_counts.png"))
    ggsave(png_file, g, width = width, height = height, dpi = dpi, bg = "white")

    if (also_pdf) {
      pdf_file <- file.path(out_dir, paste0(title, "_deg_counts.pdf"))
      ggsave(pdf_file, g, width = width, height = height, bg = "white", device = cairo_pdf)
    }
  }
  invisible(out_dir)
}


# %%
files <- list.files("../../../data/rna/pseudobulk/results/deseq2_results", pattern = "deseq2", full.names = TRUE)

batch_deg_counts(
  files,
  out_dir   = "../../../data/rna/pseudobulk/results/deg_counts_by_celltype",
  p_col     = "padj",
  es_col    = "log2FoldChange",
  term_col  = "celltype",
  facet_col = "contrast"
)
