# Author: Aishwarya Chander
# Organization: Allen Institute for Immunology

# Longitudinal frequency boxplots with subject tracking lines

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(ggplot2)
    library(ggpubr)
    library(purrr)
})

source("../utilities/functions/label-and-color-maps.R")

#' Process frequency data for plotting
#' @param df Data frame with frequency data
#' @param keep_visits Vector of visit names to keep
#' @param keep_celltypes Vector of cell types to keep
#' @param visit_col Column name for visits
#' @param celltype_col Column name for cell types
#' @return Processed data frame
process_frequency_data <- function(df,
                                   keep_visits,
                                   keep_celltypes,
                                   visit_col = "label.visitDetails",
                                   celltype_col = "aifi_plot_l3") {
    # Clean and filter data
    df <- df[df[[visit_col]] != '' & !is.na(df[[visit_col]]), ]
    
    df <- df %>%
        filter(
            !!sym(visit_col) %in% keep_visits,
            !!sym(celltype_col) %in% keep_celltypes
        ) %>%
        mutate(
            !!sym(visit_col) := factor(!!sym(visit_col), levels = keep_visits),
            celltype_label = factor(!!sym(celltype_col), levels = sort(unique(!!sym(celltype_col))))
        )
    
    return(df)
}

#' Process p-value data and calculate bracket positions
#' @param pvals_df Data frame with p-values
#' @param freq_df Frequency data frame for calculating positions
#' @param keep_visits Vector of visit names
#' @param keep_celltypes Vector of cell types
#' @param comparisons_list List of comparison pairs
#' @param value_col Column name for values (for position calculation)
#' @return Data frame with positioned p-values
process_pvalue_data <- function(pvals_df,
                                freq_df,
                                keep_visits,
                                keep_celltypes,
                                comparisons_list,
                                value_col = "cell_type_clr") {
    # Process p-values
    pvals_df <- pvals_df %>%
        mutate(
            group1 = factor(group1, levels = keep_visits),
            group2 = factor(group2, levels = keep_visits),
            celltype_label = factor(aifi_plot_l3, levels = levels(freq_df$celltype_label)),
            p.adj.signif = case_when(
                p.adj < 0.001 ~ '***',
                p.adj < 0.01 ~ '**',
                p.adj < 0.05 ~ '*',
                TRUE ~ 'ns'
            )
        ) %>%
        filter(
            p.adj.signif != 'ns',
            aifi_plot_l3 %in% keep_celltypes
        )
    
    # Filter for specified comparisons
    pvals_df <- pvals_df %>%
        filter(
            purrr::map2_lgl(group1, group2, ~
                any(purrr::map_lgl(comparisons_list, function(comp) {
                    (as.character(.x) == comp[1] & as.character(.y) == comp[2])
                }))
            )
        )
    
    # Calculate y-positions
    comparison_order <- sapply(comparisons_list, paste, collapse = '-')
    
    pvals_positioned <- pvals_df %>%
        mutate(comparison_id = paste(group1, group2, sep = '-')) %>%
        group_by(celltype_label) %>%
        mutate(
            max_y = max(freq_df[[value_col]][freq_df$celltype_label == first(celltype_label)], na.rm = TRUE),
            comparison_rank = match(comparison_id, comparison_order)
        ) %>%
        arrange(celltype_label, comparison_rank) %>%
        mutate(y.position = max_y + 0.1 + (row_number() - 1) * 0.3) %>%
        ungroup()
    
    return(pvals_positioned)
}

#' Create frequency plot with statistical annotations
#' @param df Processed frequency data
#' @param pvals_positioned Positioned p-value data
#' @param plot_params List of plot parameters
#' @return ggplot object
create_frequency_plot <- function(df,
                                  pvals_positioned = NULL,
                                  plot_params = list()) {
    # Set default parameters
    params <- list(
        n_cols = 5,
        base_font_size = 12,
        title_size = 18,
        axis_title_size = 14,
        strip_text_size = 14,
        plot_title = "Cell Type Frequencies",
        x_label = "Visit Timepoints",
        y_label = "Centered Log Ratio",
        visit_col = "label.visitDetails",
        value_col = "cell_type_clr",
        subject_col = "subject.subjectGuid"
    )
    
    # Override with user parameters
    params <- modifyList(params, plot_params)
    
    # Generate color maps from data
    visits_in_data <- levels(df[[params$visit_col]])
    tp_color_map <- get_timepoint_colors(visits_in_data)
    
    # Check if subject column exists and generate colors
    has_subjects <- params$subject_col %in% colnames(df)
    if (has_subjects) {
        subjects_in_data <- unique(df[[params$subject_col]])
        subject_color_map <- get_subject_colors(subjects_in_data)
    }
    
    # Create base plot
    p <- ggplot(df, aes(x = .data[[params$visit_col]], y = .data[[params$value_col]])) +
        geom_boxplot(
            aes(fill = .data[[params$visit_col]]),
            alpha = 0.8,
            outlier.shape = NA,
            color = 'black'
        ) +
        geom_hline(yintercept = 0, linetype = 'dashed', color = 'gray50') +
        facet_wrap(~ celltype_label, scales = 'free_y', ncol = params$n_cols) +
        labs(
            title = params$plot_title,
            x = params$x_label,
            y = params$y_label,
            fill = ''
        ) +
        theme_classic(base_size = params$base_font_size) +
        theme(
            plot.title = element_text(hjust = 0.5, face = 'bold', size = params$title_size, margin = margin(b = 10)),
            axis.title = element_text(size = params$axis_title_size, face = 'bold'),
            axis.text = element_text(size = params$base_font_size),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position = 'bottom',
            legend.title = element_text(face = 'bold', size = params$base_font_size),
            strip.text = element_text(size = params$strip_text_size, face = 'bold'),
            strip.background = element_rect(colour = NA)
        ) +
        scale_fill_manual(values = tp_color_map, guide = guide_legend(nrow = 1))
    
    # Add subject lines if available
    if (has_subjects) {
        p <- p +
            geom_line(
                aes(group = .data[[params$subject_col]], color = .data[[params$subject_col]]),
                alpha = 0.5,
                linewidth = 0.3,
                show.legend = FALSE
            ) +
            scale_color_manual(values = subject_color_map, guide = 'none')
    } else {
        p <- p + scale_color_manual(values = tp_color_map, guide = 'none')
    }
    
    # Add jittered points
    color_var <- if (has_subjects) params$subject_col else params$visit_col
    p <- p + geom_jitter(
        aes(color = .data[[color_var]]),
        width = 0.01,
        size = 1,
        stroke = 0.4,
        alpha = 0.6,
        show.legend = FALSE
    )
    
    # Add statistical annotations if provided
    if (!is.null(pvals_positioned)) {
        p <- add_statistical_brackets(p, pvals_positioned)
    }
    
    return(p)
}

#' Add statistical brackets to plot
#' @param p ggplot object
#' @param pvals_positioned Positioned p-value data
#' @param healthy_color Color for healthy comparisons
#' @return ggplot object with brackets
add_statistical_brackets <- function(p,
                                     pvals_positioned,
                                     healthy_color = "green4") {
    # Split comparisons
    healthy_comparisons <- pvals_positioned %>% filter(group2 == 'Healthy')
    other_comparisons <- pvals_positioned %>% filter(group2 != 'Healthy')
    
    # Add black brackets for non-healthy comparisons
    if (nrow(other_comparisons) > 0) {
        p <- p + stat_pvalue_manual(
            other_comparisons,
            label = 'p.adj.signif',
            tip.length = 0.01,
            size = 4,
            step.increase = 0,
            inherit.aes = FALSE
        )
    }
    
    # Add green brackets for healthy comparisons
    if (nrow(healthy_comparisons) > 0) {
        healthy_brackets <- healthy_comparisons %>%
            rowwise() %>%
            mutate(
                x1 = as.numeric(group1),
                x2 = as.numeric(group2),
                xmid = (x1 + x2) / 2
            ) %>%
            ungroup()
        
        p <- p +
            geom_segment(
                data = healthy_brackets,
                aes(x = x1, xend = x2, y = y.position, yend = y.position),
                color = healthy_color, inherit.aes = FALSE, linewidth = 0.5, alpha = 0.8
            ) +
            geom_segment(
                data = healthy_brackets,
                aes(x = x1, xend = x1, y = y.position - 0.08, yend = y.position),
                color = healthy_color, inherit.aes = FALSE, linewidth = 0.5, alpha = 0.8
            ) +
            geom_segment(
                data = healthy_brackets,
                aes(x = x2, xend = x2, y = y.position - 0.08, yend = y.position),
                color = healthy_color, inherit.aes = FALSE, linewidth = 0.5, alpha = 0.8
            ) +
            geom_text(
                data = healthy_brackets,
                aes(x = xmid, y = y.position + 0.03, label = p.adj.signif),
                color = healthy_color, inherit.aes = FALSE, size = 4, alpha = 0.8
            )
    }
    
    return(p)
}