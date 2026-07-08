# Author: Aishwarya Chander
# Organization: Allen Institute for Immunology
# Frequency lineplot with CLR-transformed values and statistical comparisons

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

source("../utilities/functions/frequency-longitudinal_boxplots.R")
source("../utilities/functions/label-and-color-maps.R")

options(warn = -1)

# Plot frequency lineplots with CLR-transformed values and significance annotations
plot_frequency_trends <- function(frequency_data_path,
                                  pvalue_data_path,
                                  keep_visits,
                                  keep_celltypes,
                                  comparisons_to_show,
                                  visit_colors = NULL,
                                  healthy_visit = "Healthy",
                                  plot_width = 15,
                                  plot_height = 6,
                                  base_size = 14,
                                  dpi = 200,
                                  n_cols = 4,
                                  title = NULL,
                                  subtitle = NULL,
                                  show_individual_points = TRUE,
                                  individual_point_alpha = 0.5,
                                  individual_point_size = 1,
                                  median_point_size = 4.5,
                                  line_color = "gray40",
                                  line_size = 1,
                                  reference_line_size = 0.8,
                                  sig_label_size = 4,
                                  sig_label_offset = 0.3,
                                  jitter_width = 0.2,
                                  show_pvals = TRUE) {
  options(repr.plot.width = plot_width, repr.plot.height = plot_height, repr.plot.res = dpi)

  # Set default visit colors if not provided
  if (is.null(visit_colors)) {
    visit_colors <- get_timepoint_colors(keep_visits)
  }

  # Load and process frequency data
  df <- fread(frequency_data_path)
  df_processed <- process_frequency_data(df, keep_visits, keep_celltypes)

  # Load and process p-values
  manual_pvals <- fread(pvalue_data_path)
  pvals_positioned <- process_pvalue_data(
    manual_pvals,
    df_processed,
    keep_visits,
    keep_celltypes,
    comparisons_to_show
  )

  # Define visits excluding healthy for trajectory
  trajectory_visits <- setdiff(keep_visits, healthy_visit)

  # Calculate summary stats per timepoint using CLR
  df_summary <- df_processed %>%
    mutate(aifi_plot_l3 = factor(aifi_plot_l3, levels = keep_celltypes)) %>%
    group_by(label.visitDetails, aifi_plot_l3) %>%
    summarise(
      median_clr = median(cell_type_clr, na.rm = TRUE),
      q25 = quantile(cell_type_clr, 0.25, na.rm = TRUE),
      q75 = quantile(cell_type_clr, 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(visit = label.visitDetails)

  # Prepare individual patient data
  df_individual <- df_processed %>%
    mutate(
      aifi_plot_l3 = factor(aifi_plot_l3, levels = keep_celltypes),
      visit = factor(label.visitDetails, levels = keep_visits)
    ) %>%
    filter(!is.na(cell_type_clr))

  # Get healthy reference values
  healthy_ref <- df_summary %>%
    filter(visit == healthy_visit) %>%
    select(aifi_plot_l3, healthy_median = median_clr)

  # Remove healthy from line data
  df_line <- df_summary %>%
    filter(visit != healthy_visit) %>%
    mutate(visit = factor(visit, levels = trajectory_visits))

  # Process stats for vs Healthy comparisons only
  stats_vs_healthy <- pvals_positioned %>%
    filter(group2 == healthy_visit) %>%
    mutate(
      visit = factor(group1, levels = trajectory_visits),
      aifi_plot_l3 = factor(aifi_plot_l3, levels = keep_celltypes),
      sig_label = case_when(
        p.adj < 0.001 ~ "***",
        p.adj < 0.01 ~ "**",
        p.adj < 0.05 ~ "*",
        TRUE ~ "n.s."
      )
    )

  # Join stats to line data for positioning
  df_plot <- df_line %>%
    left_join(healthy_ref, by = "aifi_plot_l3") %>%
    left_join(
      stats_vs_healthy %>% select(aifi_plot_l3, visit, sig_label, p.adj),
      by = c("aifi_plot_l3", "visit")
    )

  # Create base plot
  p <- ggplot(df_plot, aes(x = visit, y = median_clr))

  # Add healthy reference line
  p <- p + geom_hline(
    aes(yintercept = healthy_median),
    linetype = "dashed",
    color = visit_colors[healthy_visit],
    linewidth = reference_line_size
  )

  # Add individual patient dots if requested
  if (show_individual_points) {
    p <- p + geom_point(
      data = df_individual,
      aes(x = visit, y = cell_type_clr, color = visit),
      alpha = individual_point_alpha,
      size = individual_point_size,
      position = position_jitter(width = jitter_width, height = 0)
    )
  }

  # Add median trajectory line
  p <- p + geom_line(
    aes(group = aifi_plot_l3),
    color = line_color,
    linewidth = line_size
  )

  # Add median points
  p <- p + geom_point(
    aes(color = visit),
    size = median_point_size,
    alpha = 1
  )

  # Apply color scale
  p <- p + scale_color_manual(values = visit_colors)

  # Add significance labels if requested
  if (show_pvals) {
    p <- p + geom_text(
      data = ~ filter(.x, !is.na(sig_label)),
      aes(label = sig_label, y = q75 + sig_label_offset),
      size = sig_label_size,
      vjust = -0.5,
      color = visit_colors[healthy_visit],
      fontface = "bold"
    )
  }

  # Add faceting
  p <- p + facet_wrap(
    ~aifi_plot_l3,
    scales = "free_y",
    ncol = n_cols,
    drop = FALSE
  )

  # Add labels and theme
  p <- p + labs(
    title = title,
    subtitle = subtitle,
    x = "",
    y = "CLR-transformed frequency",
    caption = ""
  ) +
    theme_linedraw(base_size = base_size) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(face = "bold", color = "black", size = 12),
      panel.grid = element_blank(),
      # panel.grid.major.y = element_line(color = "gray90", linewidth = 0.5),
      # panel.grid.major.x = element_line(color = "gray90", linewidth = 0.5),
      legend.position = "none",
      plot.title = element_text(face = "bold")
    )

  return(p)
}
