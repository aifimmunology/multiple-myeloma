# Author: Aishwarya Chander
# Organization: Allen Institute for Immunology
# Gene expression violin plots with statistical comparisons
# ...existing code...

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(tools)
})
source("../utilities/functions/label-and-color-maps.R")
options(warn = -1)

plot_expression_violins <- function(metadata_path,
                                    files_dir,
                                    gene_list,
                                    celltypes,
                                    visits,
                                    paired_visits,
                                    unpaired_visit,
                                    visit_colors = NULL,
                                    plot_width = 10,
                                    plot_height = 10,
                                    base_size = 12,
                                    dpi = 200,
                                    title = NULL,
                                    subtitle = "",
                                    ylabs = NULL,
                                    step.increase = 0.05,
                                    tip.length = 0,
                                    bracket.nudge.y = 0.02,
                                    p.adjust.method = "BH",
                                    by_panel = TRUE,
                                    min_paired_subjects = 5,
                                    show_pvals = TRUE,
                                    show_pval_legend = TRUE,
                                    pval_colors = c("maroon", "seagreen", "gray30", "gray50"),
                                    line_size = 0.8) {
    options(repr.plot.width = plot_width, repr.plot.height = plot_height, repr.plot.res = dpi)
    if (is.null(visit_colors)) visit_colors <- get_timepoint_colors(visits)

    metadata <- fread(metadata_path, stringsAsFactors = FALSE) %>%
        filter(label.visitDetails %in% visits)

    sample_kits_needed <- unique(metadata$sample.sampleKitGuid)
    files_to_read <- file.path(files_dir, paste0(sample_kits_needed, ".csv"))

    pseudo_gex <- lapply(files_to_read, function(f) {
        if (!file.exists(f)) {
            warning(paste("File not found:", f))
            return(NULL)
        }
        df <- fread(f)
        colnames(df)[1] <- "gene"
        df %>%
            filter(gene %in% gene_list) %>%
            select(gene, any_of(celltypes)) %>%
            mutate(sample.sampleKitGuid = file_path_sans_ext(basename(f)))
    })
    pseudo_gex <- pseudo_gex[!vapply(pseudo_gex, is.null, logical(1))]

    kit_info <- metadata %>%
        select(sample.sampleKitGuid, label.visitDetails, subject.subjectGuid) %>%
        distinct()

    plot_df_long <- bind_rows(pseudo_gex) %>%
        left_join(kit_info, by = "sample.sampleKitGuid") %>%
        pivot_longer(
            cols = -c(gene, sample.sampleKitGuid, label.visitDetails, subject.subjectGuid),
            names_to = "aifi_plot_l3",
            values_to = "expression"
        ) %>%
        mutate(
            gene = factor(gene, levels = gene_list),
            label.visitDetails = factor(label.visitDetails, levels = visits)
        )

    p <- ggplot(plot_df_long, aes(x = label.visitDetails, y = expression, fill = label.visitDetails)) +
        geom_violin(scale = "width", trim = TRUE, alpha = 0.7) +
        geom_boxplot(width = 0.2, outlier.shape = 16, fill = "white", alpha = 0.8) +
        scale_fill_manual(values = visit_colors[visits]) +
        facet_grid(gene ~ aifi_plot_l3, switch = "y", scales = "free_y") +
        theme_linedraw(base_size = base_size) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            strip.text.y = element_text(angle = 0, color = "black"),
            strip.text.x = element_text(color = "black"),
            strip.background.x = element_rect(fill = "white"),
            strip.background.y = element_rect(fill = "gray90"),
            panel.grid = element_blank(),
            legend.position = "none",
            plot.title = element_text(face = "bold")
        ) +
        labs(x = "", y = ylabs, fill = "", title = title, subtitle = subtitle)

    if (show_pval_legend) {
        pval_legend_df <- data.frame(
            label.visitDetails = visits[1],
            expression = NA_real_,
            pval_group = factor(c("AdjP < 0.001", "AdjP < 0.01", "AdjP < 0.05"),
                                levels = c("AdjP < 0.001", "AdjP < 0.01", "AdjP < 0.05"))
        )

        p <- p +
            geom_point(data = pval_legend_df, aes(color = pval_group), alpha = 0, size = 0) +
            scale_color_manual(
                name = "",
                values = setNames(pval_colors[1:3], c("AdjP < 0.001", "AdjP < 0.01", "AdjP < 0.05")),
                drop = FALSE
            ) +
            guides(
                fill = "none",
                color = guide_legend(
                    override.aes = list(alpha = 1, size = 4, shape = 22, fill = pval_colors[1:3], color = pval_colors[1:3]),
                    keywidth = unit(0.8, "lines"),
                    keyheight = unit(0.8, "lines")
                )
            ) +
            theme(
                legend.position = "top",
                legend.text = element_text(size = base_size * 0.75),
                legend.key = element_rect(fill = "white", color = NA),
                legend.margin = margin(0, 0, 0, 0)
            )
    }

    if (!show_pvals) return(p)

    generate_comparisons <- function(visits, paired_visits, unpaired_visit) {
        comparisons <- list()
        paired_in_data <- intersect(paired_visits, visits)
        if (length(paired_in_data) >= 2) {
            pairs <- combn(paired_in_data, 2, simplify = FALSE)
            for (pair in pairs) {
                comparisons[[paste(pair[1], pair[2], sep = "_vs_")]] <- list(
                    group1 = pair[1], group2 = pair[2], paired = TRUE
                )
            }
        }
        if (unpaired_visit %in% visits) {
            other_visits <- setdiff(visits, unpaired_visit)
            for (v in other_visits) {
                comparisons[[paste(v, unpaired_visit, sep = "_vs_")]] <- list(
                    group1 = v, group2 = unpaired_visit, paired = FALSE
                )
            }
        }
        comparisons
    }

    comparisons <- generate_comparisons(visits, paired_visits, unpaired_visit)

    run_wilcox_test <- function(data, group1, group2, paired) {
        if (paired) {
            paired_data <- data %>%
                filter(label.visitDetails %in% c(group1, group2)) %>%
                select(subject.subjectGuid, label.visitDetails, expression) %>%
                pivot_wider(names_from = label.visitDetails, values_from = expression) %>%
                filter(!is.na(.data[[group1]]) & !is.na(.data[[group2]]))
            if (nrow(paired_data) >= min_paired_subjects) {
                return(wilcox.test(paired_data[[group1]], paired_data[[group2]], paired = TRUE)$p.value)
            }
        } else {
            vals1 <- na.omit(data$expression[data$label.visitDetails == group1])
            vals2 <- na.omit(data$expression[data$label.visitDetails == group2])
            if (length(vals1) >= 1 && length(vals2) >= 1) {
                return(wilcox.test(vals1, vals2)$p.value)
            }
        }
        NA_real_
    }

    stat_results <- plot_df_long %>%
        group_by(gene, aifi_plot_l3) %>%
        group_modify(~ {
            data <- .x
            results <- vapply(names(comparisons), function(comp_name) {
                comp <- comparisons[[comp_name]]
                run_wilcox_test(data, comp$group1, comp$group2, comp$paired)
            }, numeric(1))
            as.data.frame(t(results))
        }) %>%
        ungroup()

    adjust_pvalues <- function(stat_df, method_name, adjust_by_panel) {
        comparison_cols <- grep("^[^s].*_vs_.*$", names(stat_df), value = TRUE)
        if (adjust_by_panel) {
            stat_df %>%
                group_by(gene, aifi_plot_l3) %>%
                mutate(across(all_of(comparison_cols), ~ p.adjust(.x, method = method_name))) %>%
                ungroup()
        } else {
            for (col in comparison_cols) {
                stat_df[[col]] <- p.adjust(stat_df[[col]], method = method_name)
            }
            stat_df
        }
    }

    stat_results <- adjust_pvalues(stat_results, p.adjust.method, by_panel)

    get_bracket_color <- function(p, pval_colors) {
        case_when(
            is.na(p) ~ NA_character_,
            p < 0.001 ~ pval_colors[1],
            p < 0.01 ~ pval_colors[2],
            p < 0.05 ~ pval_colors[3],
            TRUE ~ pval_colors[4]
        )
    }

    for (comp_name in names(comparisons)) {
        stat_results[[paste0("is_significant_", comp_name)]] <-
            !is.na(stat_results[[comp_name]]) & stat_results[[comp_name]] < 0.05
        stat_results[[paste0("bracket_color_", comp_name)]] <-
            get_bracket_color(stat_results[[comp_name]], pval_colors)
    }

    y_positions <- plot_df_long %>%
        group_by(gene, aifi_plot_l3) %>%
        summarise(
            ymax = max(expression, na.rm = TRUE),
            yrange = max(expression, na.rm = TRUE) - min(expression, na.rm = TRUE),
            .groups = "drop"
        )

    stat_results <- stat_results %>%
        left_join(y_positions, by = c("gene", "aifi_plot_l3"))

    for (i in seq_along(comparisons)) {
        comp_name <- names(comparisons)[i]
        stat_results[[paste0("y_", comp_name)]] <- stat_results$ymax +
            (stat_results$yrange * bracket.nudge.y) +
            (stat_results$yrange * step.increase * (i - 1))
    }

    add_bracket <- function(p, stat_df, comp_name, group1, group2, visits) {
        sig_col <- paste0("is_significant_", comp_name)
        color_col <- paste0("bracket_color_", comp_name)
        y_col <- paste0("y_", comp_name)

        sig_data <- stat_df %>% filter(.data[[sig_col]] == TRUE)
        if (nrow(sig_data) == 0) return(p)

        visits_chr <- as.character(visits)
        x1 <- as.numeric(factor(group1, levels = visits_chr))
        x2 <- as.numeric(factor(group2, levels = visits_chr))
        if (is.na(x1) || is.na(x2)) return(p)

        y_tip_bottom <- sig_data[[y_col]] - (sig_data$yrange * tip.length)

        p +
            geom_segment(
                data = sig_data,
                aes(x = x1, xend = x2, y = .data[[y_col]], yend = .data[[y_col]]),
                inherit.aes = FALSE, size = line_size, color = sig_data[[color_col]]
            ) +
            geom_segment(
                data = sig_data,
                aes(x = x1, xend = x1, y = y_tip_bottom, yend = .data[[y_col]]),
                inherit.aes = FALSE, size = line_size, color = sig_data[[color_col]]
            ) +
            geom_segment(
                data = sig_data,
                aes(x = x2, xend = x2, y = y_tip_bottom, yend = .data[[y_col]]),
                inherit.aes = FALSE, size = line_size, color = sig_data[[color_col]]
            )
    }

    for (comp_name in names(comparisons)) {
        comp <- comparisons[[comp_name]]
        p <- add_bracket(p, stat_results, comp_name, comp$group1, comp$group2, visits)
    }

    p
}