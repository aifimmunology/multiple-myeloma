#' Sample with stratification
#'
#' This function provides a method for stratified sampling from a dataset. It allows for sampling `n` observations from each stratum defined by one or more variables. 
#' Stratification variables can include categorical variables, and NA values are treated as a separate level, ensuring that samples are representative across all levels of 
#' the stratification variables. The function is particularly useful for creating balanced samples from datasets where certain groups may be under- or over-represented.
#'
#' @param data A data.frame or tibble to sample from.
#' @param n Number of samples to draw from each stratum.
#' @param stratum Variables to stratify by, provided as a vector of column names.
#' @importFrom dplyr sample_n bind_rows
#' @importFrom stats setNames 
#' @return A data.frame or tibble containing the stratified sample.
#' @export
#' @examples
#' # Assuming `df` is your dataset and you want to sample 10 observations from each level of the `group` variable
#' df_sampled <- sampleWithStrata(df, 10, c("group"))
#' # For stratification by multiple variables, e.g., `group` and `subgroup`
#' df_sampled <- sampleWithStrata(df, 5, c("group", "subgroup"))
sampleWithStrata <- function(data, n, stratum) {
  if (!is.data.frame(data)) {
    stop("Data must be a data.frame or tibble.", call. = FALSE)
  }
  if (!all(stratum %in% names(data))) {
    stop("All stratum variables must be present in the data.", call. = FALSE)
  }
  
  # Ensure that NA values are treated as a separate level for each stratum variable
  data_with_na <- data
  for (col in stratum) {
    data_with_na[[col]] <- addNA(data[[col]])
  }
  
  # Create a unique factor that represents each combination of the stratum columns, including NA
  strata_factor <- interaction(data_with_na[stratum], drop = TRUE)
  
  # Split the data by the unique combinations of the stratum columns
  split_data <- split(data, strata_factor)
  
  # Sample n rows from each subset
  sampled_list <- lapply(split_data, function(subset) {
    if (nrow(subset) <= n) {
      return(subset)
    } else {
      return(dplyr::sample_n(subset, n))
    }
  })
  
  # Combine the sampled subsets back into one dataframe
  result <- dplyr::bind_rows(sampled_list)
  
  return(result)
}

#' Stop Execution if Condition is True
#'
#' The `stopif` function halts execution of the R script or function if the specified condition evaluates to `TRUE`.
#' This function is useful for error handling, enforcing preconditions, and ensuring that code only proceeds under
#' acceptable conditions.
#'
#' @param cond A logical condition that, when `TRUE`, causes the function to stop script execution.
#' @param message A character string that specifies the error message to be displayed upon stopping.
#' The default message is "Condition is true, stopping execution".
#'
#' @examples
#' # Stop execution if x is greater than 5
#' x <- 10
#' stopif(x > 5, "x is greater than 5, stopping execution")
#'
#' @export
stopif <- function(cond, message = "Condition is true, stopping execution") {
  if (cond) {
    stop(message)
  }
}

#' Commonest Values in Data Frame Columns
#'
#' Computes the most common values for each column in a data frame, returning the 
#' top up to 10 most common values and their frequencies, along with the class of 
#' each column. It is designed to handle numeric, factor, and character columns.
#'
#' @param df A data frame or an object convertible to a data frame, whose columns 
#' will be analyzed for the most common values.
#' @return A data frame where the first row contains the classes of the original 
#' columns and the subsequent rows contain the names and frequencies of up to 10 
#' most common values for each column, formatted as "value (frequency)".
#' @examples
#' data <- data.frame(
#'   A = sample(c("X", "Y", "Z"), 100, replace = TRUE),
#'   B = sample(1:5, 100, replace = TRUE)
#' )
#' commonest.values(data)
#' @export
commonest.values <- function(df) {
  # Error checking: Ensure df is a data frame or can be treated as one
  if (!is.data.frame(df)) {
    stop("Input must be a data frame or an object that can be coerced into one.")
  }
  
  # Compute the class for each column
  classes <- sapply(df, class)
  #if more than one class
  classes.single <- sapply(classes, paste0, collapse = "")
  # Compute the most common values and their frequencies for each column
  #common <- apply(df, 2, \(col) {
  #  most_common <- col %>% table %>% sort(decreasing = TRUE) %>% `[`(1:10) %>% paste0(names(.), " (", ., ")") #`[` pads with NA
  #})
  common <- reframe(df, across(everything(), \(col) {
    table(col) %>% 
    sort(decreasing = TRUE) %>% 
    `[`(1:10) %>% #pads with NA
    paste(names(.), " (", ., ")", sep = "")
  }))
  rbind(classes.single, as_tibble(common))
}

#' Calculate Set Operations Between Two Vectors
#'
#' This function takes two vectors and computes the union, unique elements in each vector,
#' and the intersection between them. It returns a list containing these sets.
#'
#' @param v.one A numeric or character vector.
#' @param v.two A numeric or character vector.
#' @return A list with four elements:
#'   - `all`: Union of `v.one` and `v.two`.
#'   - `unique.left`: Elements unique to `v.one`.
#'   - `unique.right`: Elements unique to `v.two`.
#'   - `common`: Intersection of `v.one` and `v.two`.
#' @examples
#' v1 <- c(1, 2, 3, 4)
#' v2 <- c(3, 4, 5, 6)
#' result <- sets(v1, v2)
#' print(result$all) # Union
#' print(result$unique.left) # Unique to v1
#' print(result$unique.right) # Unique to v2
#' print(result$common) # Intersection
#' @export
sets <- function(v.one, v.two) {
  
  # Create result list with dynamic names
  result <- list(
    all = union(v.one, v.two),
    unique.left = setdiff(v.one, v.two),
    unique.right = setdiff(v.two, v.one),
    common = intersect(v.one, v.two)
  )
  return(result)
}

#' Lengthen a Vector to a Specified Length with Fill Values
#'
#' This function extends a vector to a specified length, placing the original
#' vector's values at specified indices and filling the rest of the vector
#' with a specified fill value.
#'
#' @param vector Numeric vector to be lengthened.
#' @param indices Integer vector indicating the indices where the original
#'        values should be placed in the lengthened vector. Must be within
#'        the bounds of the specified length.
#' @param length Desired length of the output vector. Must be greater than
#'        or equal to the maximum value in `indices`.
#' @param fill Value used to fill the rest of the lengthened vector.
#'
#' @return A numeric vector of the specified length with the original values
#'         placed at the specified indices and the rest filled with the fill value.
#'
#' @examples
#' a <- c(1, 2, 3)
#' lengthen(a, indices = c(2, 4, 6), length = 10, fill = 0)
#' # Returns: c(0, 1, 0, 2, 0, 3, 0, 0, 0, 0)
#'
#' @export
lengthen <- function(vector, indices, length, fill) {
  # Error checks
  if (!is.numeric(indices) || !all(indices == floor(indices))) stop("Indices must be integer values.")
  if (length < max(indices)) stop("Length must be greater than or equal to the maximum index.")
  if (!is.numeric(length) || length != floor(length)) stop("Length must be an integer value.")
  if (length < 1) stop("Length must be at least 1.")
  if (any(indices < 1)) stop("Indices must be positive.")
  
  # Create a vector filled with the fill value
  extended_vector <- rep(fill, length)
  
  # Assign original vector values to the specified indices in the extended vector
  extended_vector[indices] <- vector
  
  return(extended_vector) 
}

#' Check Exclusivity of Values in a Data Frame Grouping
#'
#' For each group defined by a grouping column, this function checks whether the values
#' in a specified value column appear exclusively in that group and not in any other groups.
#'
#' @param df A data frame containing the data to be analyzed.
#' @param group_col A string specifying the name of the column to group by.
#' @param value_col A string specifying the name of the column containing the values to check.
#'
#' @return A named logical vector indicating for each group whether all values in the
#'         specified value column are exclusive to that group.
#' @examples
#' df <- tibble::tibble(
#'   A = c('Group1', 'Group1', 'Group2', 'Group2', 'Group3', 'Group3'),
#'   B = c(1, 2, 2, 3, 4, 5)
#' )
#' check_exclusivity(df, "A", "B")
#' 
#' @export
#' @importFrom logger log_debug
check_exclusivity <- function(df, group_col, value_col) {
  # Error checking
  if (!is.data.frame(df)) {
    stop("The input 'df' must be a data frame.")
  }
  if (!group_col %in% names(df)) {
    stop("The specified 'group_col' does not exist in the data frame.")
  }
  if (!value_col %in% names(df)) {
    stop("The specified 'value_col' does not exist in the data frame.")
  }
  
  # Splitting the data frame by the grouping variable
  groups <- split(df[[value_col]], df[[group_col]])
  
  # Checking each group against others
  exclusivity <- sapply(names(groups), function(group_name) {
    log_debug("Processing group: {group_name}")
    # The values in the current group
    current_values <- groups[[group_name]]
    
    # Values in other groups
    other_values <- unlist(groups[names(groups) != group_name])
    
    # Check if any of the current group's values are in other groups
    all(!current_values %in% other_values)
  })
  
  # Return named vector indicating exclusivity for each group
  return(exclusivity)
}

#' Load Libraries Silently with Optional Complete Suppression
#'
#' This function takes a character vector of package names and loads them into
#' the R global environment. It can suppress all messages and warnings that would
#' normally be printed to the console during the loading process, depending on
#' the value of the `suppress_all` flag. This is useful for keeping console
#' output clean, especially when loading multiple packages.
#'
#' @param packages A character vector where each element is a string
#'        representing the name of a package to be loaded.
#' @param suppress_all A logical flag indicating whether to suppress all
#'        messages and warnings (TRUE) or just messages (FALSE). Defaults to FALSE.
#'
#' @examples
#' packages_to_load <- c("dplyr", "ggplot2", "tidyr")
#' load_libraries_silently(packages_to_load, suppress_all = TRUE)
#'
#' @export
load_libraries_silently <- function(packages, suppress_all = FALSE) {
  load_function <- function(pkg) {
    if (suppress_all) {
      suppressWarnings(suppressMessages(library(pkg, character.only = TRUE)))
    } else {
      suppressMessages(library(pkg, character.only = TRUE))
    }
  }
  
  invisible(sapply(packages, load_function))
}

#' Add Metadata to SVG File
#'
#' This function adds metadata key-value pairs to an SVG file within a <metadata> element using the XML package.
#' It utilizes purrr::walk for functional iteration over the metadata elements.
#'
#' @param metadata A named list of key-value pairs to add to the SVG.
#' @param svg_filename The filepath for the SVG file to modify.
#' @return Modifies the SVG file directly by adding metadata. Does not return a value.
#' @export
#' @importFrom XML xmlParse newXMLNode xmlParent xmlRoot addChildren saveXML
#' @importFrom purrr walk
#' @examples
#' metadata_list <- list(author = "John Doe", description = "Sample SVG with metadata")
#' svg_file <- "path/to/your/svgfile.svg"
#' add_metadata_to_svg(metadata_list, svg_file)
add_metadata_to_svg <- function(metadata, svg_filename) {

  doc <- XML::xmlParse(svg_filename)
  
  metadata_node <- XML::newXMLNode("metadata")
  
  purrr::walk(names(metadata), ~{
    child_node <- XML::newXMLNode(.x, XML::newXMLTextNode(metadata[[.x]]))
    XML::addChildren(metadata_node, child_node)
  })

  root_node <- XML::xmlRoot(doc)

  XML::addChildren(root_node, metadata_node)
  
  XML::saveXML(doc, file = svg_filename)
}

#' @export
ensure_dir <- function(dir) {
    if(! dir.exists(dir)) {
        dir.create(dir, recursive=TRUE)
    }
}
