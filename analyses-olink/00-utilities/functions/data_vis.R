#' Plot a Heatmap of Column Presence Across Data Frames with Optional Transposition
#'
#' Generates a heatmap to visualize the presence (1) or absence (0) of columns across a list of data frames.
#' Each row in the heatmap represents a different column name, and each column represents a different
#' data frame from the input list. Optionally, the heatmap can be transposed to swap the axes.
#'
#' @param df_list A named list of data frames to be analyzed for column presence.
#' @param transpose A logical value indicating whether to transpose the heatmap, swapping the roles
#'        of the rows and columns. Defaults to FALSE.
#' @param ... Extra args are sent to Heatmap
#' @return An object of class `Heatmap` from the `ComplexHeatmap` package, representing the presence
#'         or absence of columns across the data frames.
#'
#' @examples
#' df1 <- data.frame(A = 1:3, B = 4:6, D = 7:9)
#' df2 <- data.frame(A = 1:3, C = 4:6, D = 7:9)
#' df3 <- data.frame(B = 1:3, C = 4:6, E = 7:9)
#' df_list <- list(df1 = df1, df2 = df2, df3 = df3)
#' plot_column_presence_heatmap(df_list, transpose = TRUE)
#'
#' @import dplyr
#' @import tidyr
#' @importFrom ComplexHeatmap Heatmap
#' @export
plot_column_presence_heatmap <- function(df_list, transpose = FALSE, ...) {
  
  # Prepare the data
  df_presence <- lapply(names(df_list), function(df_name) {
    data.frame(Column = names(df_list[[df_name]]), DF = df_name, Presence = 1)
  }) %>%
    bind_rows() %>%
    complete(DF = names(df_list), Column = unique(unlist(lapply(df_list, names)))) %>%
    replace_na(list(Presence = 0)) %>%
    pivot_wider(names_from = Column, values_from = Presence, values_fill = list(Presence = 0))

  # Convert to matrix for the Heatmap
  matrix_data <- as.matrix(select(df_presence, -DF))
  rownames(matrix_data) <- df_presence$DF

  # Transpose the matrix if requested
  if (transpose) {
    matrix_data <- t(matrix_data)
    # Swap row and column names
    #dimnames(matrix_data) <- list(colnames(matrix_data), rownames(matrix_data))
  }

  # Plot the heatmap
  Heatmap(matrix_data, name = "Presence", 
          col = c("0" = "white", "1" = "blue"),
          show_row_names = TRUE, show_column_names = TRUE, 
          cluster_rows = TRUE, cluster_columns = TRUE, 
          show_row_dend = FALSE, show_column_dend = FALSE,
          show_heatmap_legend = TRUE, ...)
}


#' Create a Presence Matrix for Heatmap Visualization
#'
#' Generates a presence matrix from a data frame using specified columns for the y-axis, x-axis,
#' and a value column, which is expected to contain a presence indicator (1 for presence).
#' Missing combinations of y-axis and x-axis are filled with 0, indicating absence.
#'
#' @param df A data frame containing the data to be visualized.
#' @param y.axis The name of the column in `df` to use as the y-axis in the heatmap.
#' @param x.axis The name of the column in `df` to use as the x-axis in the heatmap.
#' @param value The name of the column in `df` indicating presence (1) of the y-axis value in the x-axis category.
#'
#' @return A matrix suitable for heatmap visualization, with rows as y-axis values, columns as x-axis values,
#'         and cell values indicating presence (1) or absence (0).
#' @importFrom dplyr select filter mutate
#' @importFrom tidyr pivot_wider
#' @importFrom ComplexHeatmap Heatmap
#' @examples
#' df <- data.frame(
#'   Gene = rep(c("Gene1", "Gene2", "Gene3"), each = 3),
#'   Sample = rep(c("Sample1", "Sample2", "Sample3"), 3),
#'   Presence = c(1,1,0,1,0,1,1,1,1)
#' )
#' heatmap_matrix <- value_presence_matrix(df, "Gene", "Sample", "Presence")
#' ComplexHeatmap::Heatmap(heatmap_matrix)
#' 
#' @export
value_presence_matrix <- function(df, y.axis, x.axis, value, ...) {
  # Ensure the necessary packages are available

  logger::log_debug("y.axis: {y.axis}, x.axis: {x.axis}, value: {value}.")
  logger::log_debug("Column in df: {names(df)}")
  if (!all(c(y.axis, x.axis, value) %in% names(df))) {
    stop("One or more specified columns do not exist in the data frame.")
  }
  
  # Select specified columns first to handle data frames with more than three columns
  df_selected <- df %>% select(y.axis, x.axis, value)
  
  # Pivot the data frame to a wider format suitable for heatmap
  heatmap_data <- df_selected %>% 
    pivot_wider(names_from = x.axis, values_from = value, values_fill = 0)
  
  # Convert the data frame to a matrix
  heatmap_matrix <- as.matrix(select(heatmap_data, - {{y.axis}})) # Remove the first column (y-axis) before conversion
  rownames(heatmap_matrix) <- pull(heatmap_data, {{ y.axis }}) # Set row names based on y-axis
  heatmap_matrix
}

#' Save ggplot object to file with comprehensive error checking
#'
#' @param ggobject ggplot object to save.
#' @param dir Root directory under which to save the plot.
#' @param filename Name of the file without extension.
#' @param ext Extensions to save the plot as, supports "pdf", "svg", "png".
#' @param dev.opts List of device options, with 'all' for all devices or device-specific keys.
#'                 Each set of options should be a named list. Can be empty.
#' @param size List specifying width and height of the plot, defaults to 8x8.
#'
#' @examples
#' \dontrun{
#'   ggobject <- ggplot(mtcars, aes(mpg, wt)) + geom_point()
#'   plotSave(ggobject, "plots", "myplot",
#'            ext = c("pdf", "png"),
#'            dev.opts = list(all = list(units = "cm"), pdf = list(useDingbats = FALSE)),
#'            size = list(width = 10, height = 6))
#' }
#' @export
plotSave <- function(ggobject, dir, filename, ext = c("pdf", "svg", "png"),
                     dev.opts = list(), size = list(width = 8, height = 8)) {

  
  if (!is.character(dir) || !nzchar(dir)) {
    stop("dir must be a non-empty character string.")
  }
  
  if (!is.character(filename) || !nzchar(filename)) {
    stop("filename must be a non-empty character string.")
  }

  
  if (!is.list(dev.opts)) {
    stop("dev.opts must be a list.")
  }
  
  if (!is.list(size) || !all(c("width", "height") %in% names(size))) {
    stop("size must be a list containing both 'width' and 'height'.")
  }
  
  for (extension in ext) {
    save_path <- file.path(dir, extension)
    
    if (!dir.exists(save_path)) {
      dir.create(save_path, recursive = TRUE)
    }
    
    file_path <- file.path(save_path, paste0(filename, ".", extension))
    
    # Initialize options with defaults and size
    options_list <- c(list(filename = file_path, plot = ggobject, device = extension, 
                           width = size$width, height = size$height))
    
    # Merge 'all' options if any
    if (!is.null(dev.opts$all)) {
      options_list <- c(options_list, dev.opts$all)
    }
    
    # Merge device-specific options if any
    if (!is.null(dev.opts[[extension]])) {
      options_list <- c(options_list, dev.opts[[extension]])
    }
    
    do.call(ggplot2::ggsave, options_list)
  }

  ggobject
}

#' Process and Display a Table as HTML
#'
#' This function takes a frequency table, optionally sorts it, converts it to a data frame, and displays it as an HTML table.
#' Sorting can be ascending or descending but is optional.
#'
#' @param table A table object, typically generated by \code{\link[base]{table}}.
#' @param sort_column The name of the column to sort by; if NULL, no sorting is performed.
#' @param desc Boolean, whether sorting should be descending (TRUE) or ascending (FALSE). Defaults to TRUE.
#' @param var1_col_name New name to assign to the `Var1` column in the output data frame.
#' @param ... Additional arguments to pass to \code{\link[htmlTable]{htmlTable::htmlTable}}.
#'
#' @importFrom dplyr arrange desc
#' @importFrom htmlTable htmlTable
#' @importFrom IRdisplay display_html
#' @examples
#' data_vector <- sample(letters[1:4], 100, replace = TRUE)
#' tab <- table(data_vector)
#' process_and_display_table(tab, sort_column = "Freq", desc = TRUE, var1_col_name = "Category")
#' @export
process_and_display_table <- function(table, sort_column = NULL, desc = TRUE, var1_col_name = "Var1", ...) {
  # Convert table to data frame and rename the Var1 column
  df <- as.data.frame(table)
  colnames(df)[1] <- var1_col_name
  
  # Conditionally arrange the dataframe if a sort column is specified
  if (!is.null(sort_column)) {
    if (desc) {
      df <- dplyr::arrange(df, dplyr::desc(!!rlang::sym(sort_column)))
    } else {
      df <- dplyr::arrange(df, !!rlang::sym(sort_column))
    }
  }
  
  # Generate HTML table
  html_output <- htmlTable::htmlTable(df, ...)
  
  # Display HTML
  IRdisplay::display_html(html_output)
}

#' Create and Display an HTML Table from a Character Vector
#'
#' This function takes a character vector and optionally converts it into a data frame
#' with a specified column name before creating an HTML table. If no column name is
#' provided, the vector is passed directly to htmlTable. The table includes a caption
#' and can accept additional arguments for htmlTable customization.
#'
#' @param char_vector A character vector to be converted into an HTML table.
#' @param column_name Optional string specifying the column name for the data frame.
#'        If NULL, the vector is passed directly to htmlTable.
#' @param ... Additional arguments to pass to \code{\link[htmlTable]{htmlTable::htmlTable}}.
#'
#' @importFrom htmlTable htmlTable
#' @importFrom IRdisplay display_html
#' @examples
#' filenames <- c("file1.csv", "file2.csv", "data3.csv")
#' display_char_vector_as_html_table(filenames, column_name = "Filename", css.cell = "padding: 6px;")
#' display_char_vector_as_html_table(filenames, column_name = NULL, css.cell = "padding: 6px;")
#' @export
display_char_vector_as_html_table <- function(char_vector, column_name = NULL, ...) {
  # Handle the character vector based on the presence of a column name
  if (!is.null(column_name)) {
    # Convert to data frame using the specified column name and setNames for dynamic naming
    df <- setNames(data.frame(char_vector, stringsAsFactors = FALSE), column_name)
  } else {
    # Pass the vector directly to htmlTable
    df <- char_vector
  }
  
  # Create HTML table and display it
  html_output <- htmlTable::htmlTable(df, ...)
  IRdisplay::display_html(html_output)
}

#' Split Vector into N Columns and Fill Last Column with NAs
#'
#' This function takes a vector and the desired number of columns (n) as inputs.
#' It reshapes the vector into a matrix with n columns, filling in extra space in
#' the last column with `NA` values if the vector length is not divisible by n.
#'
#' @param vec A numeric or character vector that needs to be split into multiple columns.
#' @param n The number of columns into which the vector should be split.
#' @return A matrix with n columns, where the last column may contain `NA` values
#'         if the length of `vec` is not a perfect multiple of n.
#' @examples
#' vec <- 1:23
#' n_columns <- 5
#' result_matrix <- split_vector_into_columns(vec, n_columns)
#' print(result_matrix)
#' @export
split_vector_into_columns <- function(vec, n) {
  # Calculate total number of elements needed (round up to fill the matrix)
  total_elements <- ceiling(length(vec) / n) * n
  
  # Extend the vector with NAs to match the needed length
  extended_vec <- c(vec, rep(NA, total_elements - length(vec)))
  
  # Convert the extended vector to a matrix with 'n' columns
  matrix(extended_vec, nrow = total_elements / n, ncol = n, byrow = TRUE)
}

