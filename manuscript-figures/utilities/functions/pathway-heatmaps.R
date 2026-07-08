# Author: Aishwarya Chander
# Organization: Allen Institute for Immunology

# Pathway enrichment heatmaps with cell type annotations

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(stringr)
})

source("../utilities/functions/label-and-color-maps.R")

#' Select pathways by method (top_n, pattern, manual, per_cell_top)
#' @param data FGSEA results with pathway, NES columns
#' @param method Selection method
#' @param n Number of top pathways
#' @param pattern Regex pattern for filtering
#' @param pathway_list Manual pathway list
#' @param padj_threshold Significance threshold
#' @param require_significant Only include significant pathways
#' @return Character vector of pathway names
select_pathways <- function(data,
                            method = c("top_n", "pattern", "manual", "per_cell_top"),
                            n = 20,
                            pattern = NULL,
                            pathway_list = NULL,
                            padj_threshold = 0.05,
                            require_significant = FALSE) {
  # Validate method
  method <- match.arg(method)

  # Validate required columns
  if (!"pathway" %in% names(data)) stop("data must contain a 'pathway' column")
  if (!"NES" %in% names(data)) stop("data must contain a 'NES' column")
  if (!"celltype" %in% names(data) && method == "per_cell_top") {
    stop("per_cell_top requires a 'celltype' column in data")
  }

  # METHOD 1: Top N pathways by significance
  if (method == "top_n") {
    dat <- data %>% filter(!is.na(padj) & padj < padj_threshold)

    if (nrow(dat) == 0) {
      warning("No significant pathways found with padj < ", padj_threshold)
      return(character(0))
    }

    top_pathways <- dat %>%
      group_by(pathway) %>%
      summarise(
        mean_abs_nes = mean(abs(NES), na.rm = TRUE),
        min_padj = min(padj, na.rm = TRUE),
        n_sig_cells = n(),
        .groups = "drop"
      ) %>%
      arrange(desc(n_sig_cells), min_padj) %>%
      slice_head(n = n) %>%
      pull(pathway)

    return(as.character(top_pathways))

    # METHOD 2: Pattern matching
  } else if (method == "pattern") {
    if (is.null(pattern)) stop("pattern must be provided when method = 'pattern'")
    matched_pathways <- unique(data$pathway)[grepl(pattern, unique(data$pathway), ignore.case = TRUE)]

    if (length(matched_pathways) == 0) {
      warning("No pathways matched the pattern: ", pattern)
    }

    return(as.character(matched_pathways))

    # METHOD 3: Manual pathway list
  } else if (method == "manual") {
    if (is.null(pathway_list)) stop("pathway_list must be provided when method = 'manual'")

    # Case-insensitive matching while preserving order
    normalize <- function(x) toupper(trimws(x))
    data_pathways <- unique(data$pathway)
    data_norm <- normalize(data_pathways)
    list_norm <- normalize(pathway_list)

    # Preserve order from input list
    matched <- character(0)
    for (pathway in list_norm) {
      idx <- which(data_norm == pathway)
      if (length(idx) > 0) {
        matched <- c(matched, data_pathways[idx[1]])
      }
    }

    if (length(matched) == 0) {
      warning("No pathways matched from the manual list")
    }

    return(as.character(matched))

    # METHOD 4: Top N per cell type
  } else if (method == "per_cell_top") {
    dat_sel <- data

    if (require_significant) {
      dat_sel <- dat_sel %>% filter(!is.na(padj) & padj < padj_threshold)
    }

    if (nrow(dat_sel) == 0) {
      warning("No pathways available after filtering")
      return(character(0))
    }

    top_per_cell <- dat_sel %>%
      group_by(celltype, pathway) %>%
      summarise(
        mean_abs_nes = mean(abs(NES), na.rm = TRUE),
        min_padj = min(padj, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      group_by(celltype) %>%
      arrange(min_padj, desc(mean_abs_nes)) %>%
      slice_head(n = n) %>%
      ungroup() %>%
      pull(pathway) %>%
      unique()

    return(as.character(top_per_cell))
  }
}

#' Create pathway heatmap with cell type annotations
#' @param data FGSEA results with celltype, pathway, NES, padj columns
#' @param celltypes Cell types to include (NULL for all)
#' @param pathway_method Selection method: "top_n", "pattern", "manual", "per_cell_top"
#' @param n_pathways Number of top pathways
#' @param pathway_pattern Pattern for filtering
#' @param pathway_list Manual pathway list
#' @param pathway_categories Named list for faceting by category
#' @param padj_threshold Significance threshold
#' @param require_significant Require significance for selection
#' @param title Plot title
#' @param plot_width Plot width
#' @param plot_height Plot height
#' @param max_label_length Max pathway label length
#' @param base_font_size Base font size
#' @param axis_text_x_size X-axis label size
#' @param axis_text_y_size Y-axis label size
#' @param strip_text_size Facet label size
#' @param legend_title_size Legend title size
#' @param legend_text_size Legend text size
#' @param axis_text_x_angle X-axis text angle
#' @param facet_ncol Facet columns
#' @return patchwork plot object
create_pathway_heatmap <- function(data,
                                   celltypes = NULL,
                                   pathway_method = "top_n",
                                   n_pathways = 20,
                                   pathway_pattern = NULL,
                                   pathway_list = NULL,
                                   pathway_categories = NULL,
                                   padj_threshold = 0.05,
                                   require_significant = FALSE,
                                   title = NULL,
                                   plot_width = 16,
                                   plot_height = 8,
                                   max_label_length = 40,
                                   base_font_size = 16,
                                   axis_text_x_size = 12,
                                   axis_text_y_size = 12,
                                   strip_text_size = 14,
                                   legend_title_size = 10,
                                   legend_text_size = 9,
                                   axis_text_x_angle = 45,
                                   facet_ncol = NULL) {
  # Validate inputs
  required_cols <- c("celltype", "pathway", "NES", "padj")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Data missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Set celltypes if not specified
  if (is.null(celltypes)) {
    celltypes <- unique(data$celltype)
    message("Using all ", length(celltypes), " cell types found in data")
  } else {
    # Filter to cell types that exist in data
    celltypes_with_data <- unique(data$celltype)
    celltypes <- celltypes[celltypes %in% celltypes_with_data]
    message("Filtered to ", length(celltypes), " cell types that have data")
  }

  # Filter to selected celltypes
  data <- data %>% filter(celltype %in% celltypes)

  # Select pathways
  selected_pathways <- select_pathways(
    data,
    method = pathway_method,
    n = n_pathways,
    pattern = pathway_pattern,
    pathway_list = pathway_list,
    padj_threshold = padj_threshold,
    require_significant = require_significant
  )

  if (length(selected_pathways) == 0) {
    stop("No pathways selected with current criteria")
  }

  message("Selected ", length(selected_pathways), " pathways")

  # Handle pathway categories
  if (!is.null(pathway_categories)) {
    # Validate pathway_categories structure
    if (!is.list(pathway_categories) || is.null(names(pathway_categories))) {
      stop("pathway_categories must be a named list")
    }

    # Create a mapping from pathway to category, preserving order
    pathway_to_category <- character()
    category_order <- names(pathway_categories)

    for (cat_name in category_order) {
      cat_pathways <- pathway_categories[[cat_name]]
      # Match pathways case-insensitively
      normalize <- function(x) toupper(trimws(x))
      selected_norm <- normalize(selected_pathways)
      cat_norm <- normalize(cat_pathways)

      for (pathway in selected_pathways) {
        if (normalize(pathway) %in% cat_norm) {
          pathway_to_category[pathway] <- cat_name
        }
      }
    }

    # Check if all selected pathways have categories
    uncategorized <- selected_pathways[!selected_pathways %in% names(pathway_to_category)]
    if (length(uncategorized) > 0) {
      warning(
        "Some pathways don't have categories assigned: ",
        paste(head(uncategorized, 3), collapse = ", "),
        if (length(uncategorized) > 3) paste0(" (and ", length(uncategorized) - 3, " more)") else ""
      )
      # Assign them to "Other"
      for (p in uncategorized) {
        pathway_to_category[p] <- "Other"
      }
      category_order <- c(category_order, "Other")
    }

    # Filter to only pathways that are in selected_pathways
    pathway_to_category <- pathway_to_category[names(pathway_to_category) %in% selected_pathways]

    # Reorder selected_pathways based on category order and within-category input order
    ordered_pathways <- character()
    for (cat_name in category_order) {
      cat_pathways <- names(pathway_to_category)[pathway_to_category == cat_name]
      # Preserve original order within category
      if (cat_name != "Other") {
        original_order <- pathway_categories[[cat_name]]
        normalize <- function(x) toupper(trimws(x))
        cat_norm <- normalize(cat_pathways)
        orig_norm <- normalize(original_order)
        # Order by appearance in original list
        ordered_cat <- character()
        for (orig_p in original_order) {
          idx <- which(cat_norm == normalize(orig_p))
          if (length(idx) > 0) {
            ordered_cat <- c(ordered_cat, cat_pathways[idx[1]])
          }
        }
        ordered_pathways <- c(ordered_pathways, ordered_cat)
      } else {
        ordered_pathways <- c(ordered_pathways, cat_pathways)
      }
    }
    selected_pathways <- ordered_pathways
  }

  # Create complete grid for all celltype x pathway combinations
  complete_grid <- expand.grid(
    celltype = celltypes,
    pathway = selected_pathways,
    stringsAsFactors = FALSE
  )

  # Join with actual data, filling missing values
  plot_data <- complete_grid %>%
    left_join(
      data %>% filter(pathway %in% selected_pathways),
      by = c("celltype", "pathway")
    ) %>%
    mutate(
      is_na_data = is.na(NES),
      NES = ifelse(is.na(NES), 0, NES),
      padj = ifelse(is.na(padj), 1, padj),
      sig = padj < padj_threshold,
      celltype = factor(celltype, levels = celltypes),
      pathway_display = ifelse(
        nchar(pathway) > max_label_length,
        paste0(substr(pathway, 1, max_label_length), "..."),
        pathway
      ),
      l1_group = l3_to_l1[as.character(celltype)]
    ) %>%
    mutate(
      pathway_display = factor(
        pathway_display,
        levels = unique(pathway_display[order(match(pathway, selected_pathways))])
      )
    )

  # Add category information if provided
  if (!is.null(pathway_categories)) {
    plot_data <- plot_data %>%
      mutate(
        pathway_category = factor(pathway_to_category[pathway], levels = category_order),
        # Ensure pathway_display respects category order
        pathway_display = factor(pathway_display,
          levels = unique(plot_data$pathway_display[order(match(plot_data$pathway, selected_pathways))])
        )
      )
  }

  # Create L1 group annotation data
  group_data <- data.frame(
    celltype = factor(celltypes, levels = celltypes),
    l1_group = l3_to_l1[celltypes]
  )

  # Filter colors to only present groups
  present_groups <- unique(group_data$l1_group)
  group_colors <- l1_color_map[names(l1_color_map) %in% present_groups]

  # Create L1 annotation plot (side bar)
  l1_plot <- ggplot(group_data, aes(x = 1, y = celltype, fill = l1_group)) +
    geom_tile(color = NA, linewidth = 0) +
    scale_fill_manual(values = group_colors, name = "Cell Type\nGroup") +
    scale_y_discrete(drop = FALSE) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = legend_title_size, face = "bold", color = "black"),
      legend.text = element_text(size = legend_text_size, color = "black"),
      axis.text.y = element_blank(),
      plot.margin = margin(0, 0, 0, 0)
    ) +
    labs(x = NULL, y = NULL)

  # Create main heatmap
  main_plot <- ggplot(plot_data, aes(x = pathway_display, y = celltype, fill = NES)) +
    geom_tile(color = "black", linewidth = 0.2) +
    scale_fill_gradient2(
      low = "blue4",
      mid = "white",
      high = "red4",
      midpoint = 0,
      na.value = "grey90",
      name = "NES",
      limits = c(
        -max(abs(plot_data$NES), na.rm = TRUE),
        max(abs(plot_data$NES), na.rm = TRUE)
      ),
      guide = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = 1,
        barheight = 6
      )
    ) +
    scale_y_discrete(drop = FALSE) +
    theme_minimal(base_size = base_font_size) +
    theme(
      axis.text.x = element_text(
        angle = axis_text_x_angle, hjust = 1,
        size = axis_text_x_size, color = "black"
      ),
      axis.text.y = element_text(size = axis_text_y_size, color = "black"),
      axis.title = element_text(color = "black"),
      legend.position = "right",
      legend.title = element_text(size = legend_title_size, face = "bold", color = "black"),
      legend.text = element_text(size = legend_text_size, color = "black"),
      legend.box = "vertical",
      legend.box.spacing = unit(0.3, "cm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "white", color = "white"),
      strip.text = element_text(size = strip_text_size, face = "bold", color = "black"),
      panel.spacing = unit(0.2, "cm"),
      plot.title = element_blank()
    ) +
    labs(x = "", y = "")

  # Add faceting if categories are provided
  if (!is.null(pathway_categories)) {
    main_plot <- main_plot +
      facet_grid(. ~ pathway_category,
        scales = "free_x",
        space = "free_x"
      ) +
      theme(
        strip.placement = "outside",
        strip.text.x = element_text(size = strip_text_size, face = "bold", color = "black"),
        axis.text.y.left = element_text(size = axis_text_y_size, color = "black") # Only show y-axis on left
      )

    # Remove y-axis text from all but leftmost facet
    main_plot <- main_plot +
      theme(
        axis.text.y.right = element_blank() # Remove right y-axis labels
      )
  }

  # Add significance markers (white dots)
  sig_data <- plot_data %>% filter(sig == TRUE)
  if (nrow(sig_data) > 0) {
    main_plot <- main_plot +
      geom_point(
        data = sig_data,
        aes(x = pathway_display, y = celltype),
        shape = 21,
        fill = "white",
        color = "white",
        size = 2,
        stroke = 0.3,
        inherit.aes = FALSE
      )
  }

  # Combine plots with patchwork
  combined <- (main_plot + l1_plot +
    plot_layout(widths = c(0.97, 0.03), guides = "collect")) +
    plot_annotation(
      title = title,
      theme = theme(
        plot.title = if (!is.null(title)) element_text(size = 16, face = "bold", hjust = 0.5, color = "black") else element_blank(),
        legend.box = "vertical",
        legend.box.just = "left",
        legend.spacing.y = unit(0.1, "cm")
      )
    )

  # Set plot dimensions
  options(repr.plot.width = plot_width, repr.plot.height = plot_height, repr.plot.res = 200)

  return(combined)
}

#' Create faceted pathway heatmap across multiple comparisons
#' @param data_list List of FGSEA result dataframes
#' @param comparison_names Names for each comparison
#' @param celltypes Cell types to include (NULL for all)
#' @param pathway_method Selection method: "top_n", "pattern", "manual"
#' @param n_pathways Number of top pathways
#' @param pathway_pattern Pattern for filtering
#' @param pathway_list Manual pathway list
#' @param padj_threshold Significance threshold
#' @param require_significant Require significance for selection
#' @param title Plot title
#' @param plot_width Plot width
#' @param plot_height Plot height
#' @param max_label_length Max pathway label length
#' @param facet_ncol Facet columns
#' @param base_font_size Base font size
#' @param axis_text_x_size X-axis label size
#' @param axis_text_y_size Y-axis label size
#' @param strip_text_size Facet label size
#' @param legend_title_size Legend title size
#' @param legend_text_size Legend text size
#' @param axis_text_x_angle X-axis text angle
#' @return patchwork plot object
create_faceted_pathway_heatmap <- function(data_list,
                                           comparison_names,
                                           celltypes = NULL,
                                           pathway_method = "top_n",
                                           n_pathways = 20,
                                           pathway_pattern = NULL,
                                           pathway_list = NULL,
                                           padj_threshold = 0.05,
                                           require_significant = FALSE,
                                           title = NULL,
                                           plot_width = 18,
                                           plot_height = 8,
                                           max_label_length = 30,
                                           facet_ncol = 4,
                                           base_font_size = 16,
                                           axis_text_x_size = 12,
                                           axis_text_y_size = 12,
                                           strip_text_size = 11,
                                           legend_title_size = 10,
                                           legend_text_size = 9,
                                           axis_text_x_angle = 45) {
  # Validate inputs
  if (length(data_list) != length(comparison_names)) {
    stop("data_list and comparison_names must have the same length")
  }

  if (length(data_list) == 0) {
    stop("data_list cannot be empty")
  }

  # Get all unique celltypes if not specified
  if (is.null(celltypes)) {
    celltypes <- unique(unlist(lapply(data_list, function(x) unique(x$celltype))))
    message("Using all ", length(celltypes), " cell types found across datasets")
  } else {
    # Filter to cell types that exist in data
    celltypes_with_data <- unique(unlist(lapply(data_list, function(x) unique(x$celltype))))
    celltypes <- celltypes[celltypes %in% celltypes_with_data]
    message("Filtered to ", length(celltypes), " cell types that have data across datasets")
  }

    celltype_order <- names(l3_to_l1)
    celltypes <- celltypes[order(match(celltypes, celltype_order))]
    message("Reordered cell types based on L3 to L1 hierarchy")
                                                
  # Combine all data for pathway selection (ensures consistency across facets)
  combined_for_selection <- do.call(rbind, lapply(seq_along(data_list), function(i) {
    data_list[[i]] %>%
      filter(celltype %in% celltypes) %>%
      mutate(comparison = comparison_names[i])
  }))

  # Select pathways from combined data
  selected_pathways <- select_pathways(
    combined_for_selection,
    method = pathway_method,
    n = n_pathways,
    pattern = pathway_pattern,
    pathway_list = pathway_list,
    padj_threshold = padj_threshold,
    require_significant = require_significant
  )

  if (length(selected_pathways) == 0) {
    stop("No pathways selected with current criteria")
  }

  message("Selected ", length(selected_pathways), " pathways from combined data")

  # Prepare data for each comparison with complete grid
  prepared_data_list <- lapply(seq_along(data_list), function(i) {
    data <- data_list[[i]] %>%
      filter(celltype %in% celltypes)

    # Create complete grid (all celltype x pathway combinations)
    complete_grid <- expand.grid(
      celltype = celltypes,
      pathway = selected_pathways,
      stringsAsFactors = FALSE
    )

    # Join with actual data, filling missing values
    plot_data <- complete_grid %>%
      left_join(
        data %>% filter(pathway %in% selected_pathways),
        by = c("celltype", "pathway")
      ) %>%
      mutate(
        is_na_data = is.na(NES),
        NES = ifelse(is.na(NES), 0, NES),
        padj = ifelse(is.na(padj), 1, padj),
        sig = padj < padj_threshold,
        comparison = comparison_names[i],
        celltype = factor(celltype, levels = celltypes),
        pathway_display = ifelse(
          nchar(pathway) > max_label_length,
          paste0(substr(pathway, 1, max_label_length), "..."),
          pathway
        ),
        l1_group = l3_to_l1[as.character(celltype)]
      )

    return(plot_data)
  })

  # Combine all prepared data
  combined_plot_data <- do.call(rbind, prepared_data_list) %>%
    mutate(comparison = factor(comparison, levels = comparison_names))

  # Create L1 group annotation data
  group_data <- data.frame(
    celltype = factor(celltypes, levels = celltypes),
    l1_group = l3_to_l1[celltypes]
  )

  # Filter colors to only present groups
  present_groups <- unique(group_data$l1_group)
  group_colors <- l1_color_map[names(l1_color_map) %in% present_groups]

  # Create L1 annotation plot (side bar)
  l1_plot <- ggplot(group_data, aes(x = 1, y = celltype, fill = l1_group)) +
    geom_tile(color = NA, linewidth = 0) +
    scale_fill_manual(values = group_colors, name = "Cell Type\nGroup") +
    scale_y_discrete(drop = FALSE) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = legend_title_size, face = "bold", color = "black"),
      legend.text = element_text(size = legend_text_size, color = "black"),
      axis.text.y = element_blank(),
      plot.margin = margin(0, 0, 0, 0)
    ) +
    labs(x = NULL, y = NULL)

  # Create main faceted heatmap
  main_plot <- ggplot(combined_plot_data, aes(x = pathway_display, y = celltype, fill = NES)) +
    geom_tile(color = "black", linewidth = 0.2) +
    scale_fill_gradient2(
      low = "blue4",
      mid = "white",
      high = "red4",
      midpoint = 0,
      na.value = "grey90",
      name = "NES",
      limits = c(
        -max(abs(combined_plot_data$NES), na.rm = TRUE),
        max(abs(combined_plot_data$NES), na.rm = TRUE)
      ),
      guide = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = 1,
        barheight = 6
      )
    ) +
    scale_y_discrete(drop = FALSE) +
    facet_wrap(~comparison, ncol = facet_ncol, scales = "free_x") +
    theme_minimal(base_size = base_font_size) +
    theme(
      axis.text.x = element_text(
        angle = axis_text_x_angle, hjust = 1,
        size = axis_text_x_size, color = "black"
      ),
      axis.text.y = element_text(size = axis_text_y_size, color = "black"),
      axis.title = element_text(color = "black"),
      legend.position = "right",
      legend.title = element_text(size = legend_title_size, face = "bold", color = "black"),
      legend.text = element_text(size = legend_text_size, color = "black"),
      legend.box = "vertical",
      legend.box.spacing = unit(0.3, "cm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "white", color = "white"),
      strip.text = element_text(size = strip_text_size, face = "bold", color = "black"),
      panel.spacing = unit(0.2, "cm"),
      plot.margin = margin(5, 5, 5, 5),
      plot.title = element_blank()
    ) +
    labs(x = "", y = "")

  # Add significance markers (white dots)
  sig_data <- combined_plot_data %>% filter(sig == TRUE)
  if (nrow(sig_data) > 0) {
    main_plot <- main_plot +
      geom_point(
        data = sig_data,
        aes(x = pathway_display, y = celltype),
        shape = 21,
        fill = "white",
        color = "white",
        size = 2,
        stroke = 0.3,
        inherit.aes = FALSE
      )
  }

  # Combine plots with patchwork
  combined <- (main_plot + l1_plot +
    plot_layout(widths = c(0.97, 0.03), guides = "collect")) +
    plot_annotation(
      title = title,
      theme = theme(
        plot.title = if (!is.null(title)) element_text(size = 16, face = "bold", hjust = 0.5, color = "black") else element_blank(),
        legend.box = "vertical",
        legend.box.just = "left",
        legend.spacing.y = unit(0.1, "cm")
      )
    )

  # Set plot dimensions
  options(repr.plot.width = plot_width, repr.plot.height = plot_height, repr.plot.res = 200)

  return(combined)
}