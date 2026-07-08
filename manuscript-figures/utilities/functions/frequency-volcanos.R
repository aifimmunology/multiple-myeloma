# Author: Aishwarya Chander
# Organization: Allen Institute for Immunology

# Volcano plots for cell type frequency changes

suppressPackageStartupMessages({
    library(tidyverse)
    library(ggrepel)
    library(data.table)
})

source("../utilities/functions/label-and-color-maps.R")

#' Unpaired Wilcoxon test for frequency comparison
#' @param freq_df Frequency data frame
#' @param g1 First group name
#' @param g2 Second group name
#' @param freq_data Column name for frequency values
#' @param label_col Column name for cell type labels
#' @return Data frame with test results
test_wilcox_unpaired <- function(freq_df, g1, g2,
                                 freq_data = 'cell_type_clr',
                                 label_col = 'aifi_plot_l3') {
  freq_df %>%
    filter(label.visitDetails %in% c(g1, g2)) %>%
    distinct(subject.subjectGuid, .data[[label_col]], label.visitDetails, .data[[freq_data]]) %>%
    group_by(.data[[label_col]]) %>%
    reframe({
      x <- .data[[freq_data]][label.visitDetails == g1]
      y <- .data[[freq_data]][label.visitDetails == g2]
      wt <- wilcox.test(x, y, paired = FALSE, exact = TRUE)
      tibble(
        n1 = length(x),
        n2 = length(y),
        statistic = unname(wt$statistic),
        p = wt$p.value
      )
    }) %>%
    mutate(
      group1 = g1,
      group2 = g2,
      .before = 1
    ) %>%
    arrange(.data[[label_col]]) %>%
    mutate(p.adj = p.adjust(p, method = 'fdr'))
}

#' Calculate median difference effect size
#' @param freq_df Frequency data frame
#' @param g1 First group name
#' @param g2 Second group name
#' @param freq_data Column name for frequency values
#' @param label_col Column name for cell type labels
#' @return Data frame with effect sizes
calc_effect_size <- function(freq_df, g1, g2, freq_data, label_col) {
  freq_df %>%
    filter(label.visitDetails %in% c(g1, g2)) %>%
    distinct(subject.subjectGuid, .data[[label_col]], label.visitDetails, .data[[freq_data]]) %>%
    group_by(.data[[label_col]]) %>%
    summarise(
      diff = median(.data[[freq_data]][label.visitDetails == g2], na.rm = TRUE) -
             median(.data[[freq_data]][label.visitDetails == g1], na.rm = TRUE),
      .groups = 'drop'
    )
}

#' Create frequency change volcano plot
#' @param freq_df Frequency data frame
#' @param deg_res Wilcoxon test results from test_wilcox_unpaired
#' @param freq_data Column name for frequency values
#' @param label_col Column name for cell type labels
#' @param alpha Significance threshold
#' @param effect_size Minimum effect size threshold
#' @param title Plot title
#' @return ggplot object
plot_freq_volcano <- function(freq_df, deg_res,
                freq_data = "cell_type_clr",
                label_col = "aifi_plot_l3",
                alpha = 0.05,
                effect_size = 0,
                title = NULL) {

  # Extract group names
  g1 <- unique(deg_res$group1)[1]
  g2 <- unique(deg_res$group2)[1]

  # Get colors from centralized map
  col_g1 <- get_timepoint_colors(g1)
  col_g2 <- get_timepoint_colors(g2)
  
  # Calculate effect sizes and prepare plotting data
  eff <- calc_effect_size(freq_df, g1, g2, freq_data, label_col)
  
  res_sig <- deg_res %>%
  left_join(eff, by = label_col) %>%
  mutate(
    neglog = -log10(pmax(p.adj, .Machine$double.eps)),
    sig = p.adj <= alpha & abs(diff) > effect_size,
    dir = case_when(
    sig & diff > 0 ~ paste("Up in", g2),
    sig & diff < 0 ~ paste("Up in", g1),
    TRUE ~ "Not sig"
    ),
    dir = factor(dir, levels = c(paste("Up in", g1), "Not sig", paste("Up in", g2))),
    label_wrapped = str_wrap(.data[[label_col]], width = 25)
  )
  
  # Define color mapping
  col_map <- c(
  setNames(col_g1, paste("Up in", g1)),
  setNames("grey70", "Not sig"),
  setNames(col_g2, paste("Up in", g2))
  )
  
  # Create plot
  p <- ggplot(res_sig, aes(x = diff, y = neglog, color = dir)) +
  geom_point(size = 3, alpha = 0.6) +
  geom_text_repel(
    data = filter(res_sig, sig),
    aes(label = label_wrapped),
    size = 6.5,
    max.overlaps = 4,
    min.segment.length = 0,
    show.legend = FALSE
  ) +
  geom_hline(yintercept = -log10(alpha), linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = c(-effect_size, effect_size), linetype = "dashed", alpha = 0.5) +
  scale_color_manual(values = col_map, drop = FALSE) +
  labs(
    title = title %||% paste0(g1, " vs ", g2),
    x = "Median difference",
    y = "-log10(FDR)",
    color = NULL
  ) +
  theme_classic(base_size = 18) +
  theme(
    plot.title = element_text(size = 17, face = "bold", margin = margin(b = 5)),
    legend.position = "bottom",
    legend.text = element_text(margin = margin(r = 10, l = 0))
  )
 
  return(p)
}