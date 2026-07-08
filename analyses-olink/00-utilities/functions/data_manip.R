#' Map Column Names and Remove Columns in a Data Frame
#'
#' This function applies a mapping to rename and/or remove columns in a given data frame. 
#' Columns can be removed by mapping their names to `NULL` and renamed by mapping their names 
#' to new character values. The function checks for input validity, handles potential naming 
#' conflicts, and logs actions for transparency.
#'
#' @param df A data frame object on which column mapping is to be applied.
#' @param mapping A named list where the names are the original column names in `df` and the 
#' values are the new names for those columns. A value of `NULL` indicates that the column 
#' should be removed.
#' @param delete.missng should names that are in df but not in mapping by default be deleted?
#'
#' @return A modified data frame with columns renamed and/or removed according to `mapping`.
#' @export
#'
#' @examples
#' df <- data.frame(a = 1:3, b = 4:6, c = 7:9)
#' mapping <- list(a = "alpha", b = NULL, c = "gamma")
#' new_df <- column_map(df, mapping)
#' # new_df will have columns "alpha" and "gamma", with "b" removed
#'
#' @importFrom dplyr select rename_with
#' @importFrom logger log_info log_debug
column_map <- function(df, mapping, delete.missing = F) {
  # Check if the input is a data frame
  if (!inherits(df, "data.frame")) {
    stop("df must inherit class data.frame", call. = FALSE)
  }

  stopifnot(length(mapping) > 0)
  
  if(delete.missing) {
    missing.in.mapping <- setdiff(names(df), names(mapping))
    if(length(missing.in.mapping) > 0) {
      logger::log_info("Added column for deletion due to delete.missing = T: {missing.in.mapping}")
      missing.to.delete <- setNames(vector("list", length(missing.in.mapping)), missing.in.mapping)
      missing.to.delete <- lapply(missing.to.delete, function(x) NULL)
      mapping <- c(mapping, missing.to.delete)
      }
  }
  
  #columns indicated for deletion in mapping
  mapping_delete_cols <- names(mapping)[sapply(mapping, is.null)]

  #columns indicated for renaming in the mapping 
  mapping_rename_cols <- setdiff(names(mapping), mapping_delete_cols)

  mapping_rename_type_correct <- sapply(mapping[mapping_rename_cols], is.character)
  if(! all(mapping_rename_type_correct)) {
      logger::log_fatal("Found invalid type in mapping amoong columns for renaming cause by {mapping_rename_cols[! mapping_rename_type_correct]}")
      
      stop("The mapping can only contain NULL to indicate delete or character to indicate new name")
  }
    

  #columns to be deleted in df
  df_delete_cols <- names(df)[names(df) %in% mapping_delete_cols]
  logger::log_info("Will delete column: {df_delete_cols}")
    
  # Apply renaming and removing in a tidy way
  df_cols_deleted <- df %>%
    dplyr::select(-all_of(df_delete_cols)) 

    
    
  #columns that exist in df to rename
  df_rename_cols <- mapping_rename_cols[mapping_rename_cols %in% names(df)]
    
  # Create a named vector for setNames function
  rename_vector <- unlist(mapping[df_rename_cols])

  # Remove items where there is no renaming, i.e. old an new name is the same
  rename_vector_action <- rename_vector[names(rename_vector) != rename_vector]

  # If the new name of a column duplicates another, e.g. if both NPX and NPX_bridged -> NPX exist, duplicates will arise. 
  # delete the old column before renaming to avoid that
  rename_conflicts <- rename_vector_action %in% names(df_cols_deleted)
  if(any(rename_conflicts)) {
    cols_conflict <- rename_vector_action[rename_conflicts]
    logger::log_info("Rename Conflict: Removing old column {cols_conflict} before renaming {names(rename_vector_action)[rename_conflicts]}")
    df_cols_deleted <- select(df_cols_deleted, -all_of(cols_conflict))
  }
  
  logger::log_debug("Will rename column {names(rename_vector_action)} to {rename_vector_action}")

  df_mapped <- df_cols_deleted %>% dplyr::rename_with(~rename_vector_action[.x], .cols = all_of(names(rename_vector_action)))
  
  return(df_mapped)
}

#' Check One-to-One Mapping Between Columns
#'
#' This function examines whether there is a one-to-one mapping between values in two specified columns
#' of a data frame. If a value in the reference column maps to more than one value in the target column,
#' those mappings are returned in a long dataframe for inspection.
#'
#' @param df A data frame to be inspected.
#' @param reference_column The name of the reference column as a string.
#' @param target_column The name of the target column as a string.
#'
#' @return A data frame listing reference values that map to more than one target value,
#'         along with their corresponding target values for inspection. Returns `NULL` if all
#'         reference values map to exactly one target value.
#'
#' @examples
#' df <- data.frame(
#'   ref_col = c("A", "A", "B", "B", "C"),
#'   target_col = c(1, 1, 2, 3, 4)
#' )
#' check_one_to_one_mapping(df, "ref_col", "target_col")
#'
#' @importFrom dplyr select distinct group_by summarise filter n_distinct pull
#' @importFrom tidyr unite
#' @export
check_one_to_one_mapping <- function(df, reference_column, target_column, .id=NULL) {
  # Ensure df is a data frame
  if (!inherits(df, "data.frame")) {
    stop("The input `df` must be a data frame.", call. = FALSE)
  }
  
  # Check if the specified columns and .id (if provided) exist in the data frame
  required_columns <- c(reference_column, target_column)
  if (!is.null(.id)) {
    required_columns <- c(required_columns, .id)
    if (!(.id %in% names(df))) {
      stop("The .id column must exist in the data frame.", call. = FALSE)
    }
  }
  
  if (!all(required_columns %in% names(df))) {
    stop("Both reference and target columns (and .id if specified) must exist in the data frame.", call. = FALSE)
  }
  
  # Proceed with checking one-to-one mapping
  df_selected <- df %>%
    select(all_of(required_columns)) %>%
    distinct()
  
  # Count the number of unique target values for each reference value
  mapping_counts <- df_selected %>%
    group_by(.data[[reference_column]]) %>%
    summarise(n_target_values = n_distinct(.data[[target_column]]), .groups = 'drop')
  
  # Identify reference values with more than one corresponding target value
  multiple_mappings <- mapping_counts %>%
    filter(n_target_values > 1) %>%
    pull(.data[[reference_column]])
  
  if (length(multiple_mappings) == 0) {
    logger::log_info("All values in the reference column map to exactly one value in the target column.")
    return(NULL)
  } else {
    message <- paste("The following values in the", reference_column, "column have multiple mappings in the", target_column, "column:", toString(multiple_mappings))
    logger::log_info(message)
    
    # Return a dataframe for inspection of mappings, including .id if specified
    df_for_inspection <- df_selected %>%
      filter(.data[[reference_column]] %in% multiple_mappings) %>%
      arrange(match(.data[[reference_column]], multiple_mappings))
    
    return(df_for_inspection)
  }
}

#' Test for One-to-One Mapping Between Two Columns
#'
#' This function checks if there is a one-to-one mapping between two specified columns in a dataframe,
#' returning a boolean value. It utilizes `check_one_to_one_mapping` to perform the check.
#' `TRUE` is returned if a one-to-one mapping exists; otherwise, `FALSE` is returned.
#'
#' For details on the parameters `df`, `reference_column`, `target_column`, and `.id`, 
#' see the documentation for `check_one_to_one_mapping`.
#'
#' @return A boolean value; `TRUE` if a one-to-one mapping exists between the `reference_column` and `target_column`, `FALSE` otherwise.
#'
#' @seealso \code{\link{check_one_to_one_mapping}} for details on the parameters and the underlying logic used to perform the mapping check.
#'
#' @examples
#' df <- data.frame(
#'   key = c("a", "b", "c", "a"),
#'   value = c(1, 2, 3, 1)
#' )
#' has_one_to_one_mapping(df, "key", "value")
#' # Expected output: FALSE, because "a" maps to 1 more than once
#'
#' @export
has_one_to_one_mapping <- function(df, reference_column, target_column, .id=NULL) {
  res <- check_one_to_one_mapping(df=df, reference_column=reference_column, target_column=target_column, .id=.id)
  is.null(res) # NULL indicates one-to-one mapping
}

                                  
#' Get Common Columns Across Data Frames
#'
#' Identifies and returns the common column names present in all data frames within a given list of data frames.
#' 
#' @param df_list A list of data frames for which to find common columns. Each element of the list must be a data frame.
#' 
#' @return A character vector containing the names of columns that are present in all data frames in the input list.
#' 
#' @examples
#' df1 <- data.frame(a = 1, b = 2, c = 3)
#' df2 <- data.frame(a = 4, b = 5, d = 6)
#' df_list <- list(df1, df2)
#' common_cols <- get_common_columns(df_list)
#'
#' @export
#'
#' @importFrom purrr map reduce
#' @importFrom dplyr select
get_common_columns <- function(df_list) {
  # Validate input: ensure each element is a data frame
  walk(df_list, \(df) {
    if (!inherits(df, "data.frame")) {
      stop("Each element of `df_list` must be a data frame.", call. = FALSE)
    }
  })
  
  # Identify common columns
  common_columns <- df_list %>%
    map(names) %>%
    reduce(intersect)
  
  return(common_columns)
}

#' Stack Data Frames by Common Columns
#'
#' Stacks a list of data frames by common columns, optionally limiting to a subset of columns.
#' Logs information about dropped columns if there are any columns not common across all data frames.
#' Ensures that all input data frames are aligned by common column names before stacking.
#'
#' @param df_list A list of data frames to be stacked. Each element of the list must be a data frame.
#' @param subset An optional character vector specifying a subset of common columns to retain.
#'               If NULL, all common columns are retained.
#' @param id Passed ont .id in bind_rows
#' @return A data frame resulting from binding the rows of all input data frames, limited to common columns
#'         (and further limited to a specified subset of columns if provided).
#'
#' @examples
#' df1 <- data.frame(a = 1, b = 2, c = 3)
#' df2 <- data.frame(a = 4, b = 5, d = 6)
#' df3 <- data.frame(a = 7, b = 8, c = 9)
#' df_list <- list(df1, df2, df3)
#' stacked_df <- stack_dfs(df_list)
#'
#' @export
#'
#' @importFrom purrr map reduce
#' @importFrom dplyr select all_of bind_rows
#' @importFrom logger log_info
stack_dfs <- function(df_list, id=NULL) {
  walk(df_list, \(df) {
    if (!inherits(df, "data.frame")) {
      stop("Each element of `df_list` must be a data frame.", call. = FALSE)
    }
  })
  
  common_columns <- get_common_columns(df_list)
  logger::log_info("Found common column: {common_columns}")
  all_columns <- reduce(df_list %>% map(names), union)
  disjoint_columns <- setdiff(all_columns, common_columns)
    
  if (length(disjoint_columns) > 0) {
    logger::log_info(paste("Dropping column(s):", paste(disjoint_columns, collapse = ", ")))
  }

  
  stopifnot(length(common_columns) > 0)
  
  df_list_common <- map(df_list, ~select(.x, all_of(common_columns)))

#  if (!is.null(subset)) {
#    stopifnot(all(subset %in% common_columns))
#    df_list_common <- map(df_list_common, ~select(.x, all_of(subset)))
#  }
  
  if(!is.null(id)) print(id)
  dplyr::bind_rows(df_list_common, .id = id)
}

#' Justify Substring to Right of Delimiter
#'
#' This function takes a character vector and a delimiter, finds the maximum length from the start
#' of each string to the delimiter (plus one for spacing), and then pads each string with whitespace
#' to justify the substring to the right of the delimiter into a column.
#'
#' @param strings A character vector to be processed.
#' @param delimiter A character string representing the delimiter used to split the strings.
#'
#' @return A character vector with each string padded with whitespace to align substrings
#' to the right of the delimiter.
#'
#' @examples
#' strings <- c("key:value", "longerkey:value2", "short:value3")
#' justify_right_of_delimiter(strings, ":")
#'
#' @export
justify_right_of_delimiter <- function(strings, delimiter) {
  # Find the position of the delimiter in each string
  delimiter_positions <- sapply(strings, function(s) {
    match_position <- regexpr(delimiter, s)
    if (match_position < 0) return(nchar(s) + 1) # No delimiter found, use string length
    return(match_position + 1) # Adjust for zero-index
  })
  
  # Determine the maximum length to the delimiter across all strings
  max_len <- max(delimiter_positions)
  
  # Pad each string with whitespace to justify the substring
  justified_strings <- sapply(strings, function(s) {
    current_len <- regexpr(delimiter, s) + 1
    if (current_len < 0) current_len <- nchar(s) + 1
    padding_needed <- max_len - current_len
    #paste0(paste(rep(" ", padding_needed), collapse = ""), s)
    gsub(delimiter, paste0(rep(" ", padding_needed), collapse=""), s)
      
  })

  return(justified_strings)
}


#' Identify Varying Columns Among Duplicates
#'
#' This function takes a dataframe that contains duplicate rows based on a specified column
#' and identifies columns that have varying data among these duplicates. It is particularly
#' useful for analyzing the result of `janitor::get_dupes()` to find out which columns contribute
#' to the distinction among duplicated entries.
#'
#' @param df A dataframe with duplicates, ideally obtained from `janitor::get_dupes()`.
#' @param dupe_column The name of the column that contains the duplicate identifiers as a string.
#'
#' @return A dataframe that includes the original columns specified plus a new column `cols_to_blame`
#'         indicating which columns vary among the duplicates for each unique value in `dupe_column`.
#'         It also performs a join to include only those rows that have been identified as duplicates
#'         and adds information about varying columns.
#'
#' @examples
#' # Assuming `your_dataframe` is a dataframe obtained from `janitor::get_dupes()`
#' # and you want to find varying columns for duplicates in "sample.sampleKitGuid":
#' result <- blame.dupes(your_dataframe, "sample.sampleKitGuid")
#'
#' @importFrom dplyr group_by summarise across filter mutate select distinct left_join
#' @importFrom tidyr pivot_longer
#' @importFrom rlang sym
#' @export
#' @rdname blame.dupes                                
blame.dupes <- function(df, dupe_column) {
  # Ensure the column is treated as a symbol for tidy evaluation
  dupe_col_sym <- rlang::sym(dupe_column)
  
  # Find varying columns among duplicates
  blame.df <- df %>%
    group_by(!!dupe_col_sym) %>%
    summarise(across(everything(), ~n_distinct(.) > 1), .groups = "drop") %>%
    pivot_longer(-!!dupe_col_sym, names_to = "column", values_to = "blame") %>%
    filter(blame) %>%
    group_by(!!dupe_col_sym) %>%
    mutate(cols_to_blame = paste0(column, collapse = ", ")) %>%
    distinct()
  
  # Select those columns that vary among duplicates
  descriptors.dupes.causative <- df %>%
    select(!!dupe_col_sym, dupe_count, all_of(unique(blame.df$column)))
  
  # Join back and include a column to indicate which columns are varying
  descriptors.dupes.with.causative.cols <- blame.df %>%
    select(!!dupe_col_sym, cols_to_blame) %>%
    distinct() %>%
    left_join(descriptors.dupes.causative, by = dupe_column)
  
  return(descriptors.dupes.with.causative.cols)
}

#' Equalize Lengths of Vectors in a List
#'
#' @param vectors A list of vectors to be equalized in length.
#' @param fill Value used to fill the extended parts of the vectors. Default is NA.
#' @return A list of vectors, all of equal length, filled with the specified value if extended.
#' @examples
#' list_of_vectors <- list(c(1, 2, 3), c(4, 5), c(6, 7, 8, 9))
#' equalized_vectors <- equalize_vector_lengths(list_of_vectors)
equalize_vector_lengths <- function(vectors, fill = NA) {
  if (!is.list(vectors)) stop("The 'vectors' argument must be a list.")
  if (!all(sapply(vectors, is.vector))) stop("All elements of 'vectors' must be vectors.")
  
  max_length <- max(sapply(vectors, length))
  lapply(vectors, function(x) {
    length(x) <- max_length
    x[is.na(x)] <- fill
    return(x)
  })
}

#' Convert a List of Vectors to a Data Frame with Equalized Lengths
#'
#' @param list_of_vectors A list of vectors of potentially varying lengths.
#' @param fill Value used to fill the extended parts of the vectors if they are extended to match the longest vector. Default is NA.
#' @return A data frame where each vector from the list becomes a column.
#' @export
#' @examples
#' list_of_vectors <- list(c(1, 2, 3), c(4, 5), c(6, 7, 8, 9))
#' df <- list2df(list_of_vectors)
list2df <- function(list_of_vectors, fill = NA) {
  if (!is.list(list_of_vectors)) stop("The input must be a list of vectors.")
    
  equalized_vectors <- equalize_vector_lengths(list_of_vectors, fill = fill)
  df <- tibble::as_tibble(equalized_vectors)

   
  return(df)
}

#' Equalize Lengths of Data Frames in a List
#'
#' Pads data frames in the list to make them all have the same number of rows, matching the longest data frame.
#' Rows added to shorter data frames are filled with NA or a specified value.
#'
#' @param list_of_dfs A list of data frames to be equalized in row number.
#' @param fill Value used to fill the extended parts of the data frames. Default is NA.
#' @return A list of data frames, all of equal row number, filled with the specified value if extended.
#' @examples
#' df1 <- data.frame(x = 1:3, y = letters[1:3])
#' df2 <- data.frame(x = 1:5, y = letters[1:5])
#' equalized_dfs <- equalize_df_lengths(list(df1, df2))
equalize_df_lengths <- function(list_of_dfs, fill = NA) {
  if (!is.list(list_of_dfs)) stop("The input must be a list.")
  if (!all(sapply(list_of_dfs, is.data.frame))) stop("All elements of the list must be data frames.")
  
  max_rows <- max(sapply(list_of_dfs, nrow))
  lapply(list_of_dfs, function(df) {
    n <- nrow(df)
    if (n < max_rows) {
      extra_rows <- max_rows - n
      df <- rbind(df, as.data.frame(matrix(rep(fill, extra_rows * ncol(df)), nrow = extra_rows, 
                                           dimnames = list(NULL, names(df)))))
    }
    return(df)
  })
}

#' Convert a List of Data Frames to a Single Data Frame with Equalized Row Numbers
#'
#' @param list_of_dfs A list of data frames of potentially varying row numbers.
#' @param fill Value used to fill the rows added to data frames to match the number of rows of the longest data frame. Default is NA.
#' @return A single data frame by binding the columns of all data frames in the list, after equalizing their row numbers.
#' @export
#' @examples
#' df1 <- data.frame(x = 1:3, y = letters[1:3])
#' df2 <- data.frame(x = 1:5, y = letters[1:5])
#' combined_df <- list_dfs2df(list(df1, df2))
list_dfs2df <- function(list_of_dfs, fill = NA) {
  if (!is.list(list_of_dfs)) stop("The input must be a list of data frames.")
  
  equalized_dfs <- equalize_df_lengths(list_of_dfs, fill = fill)
  combined_df <- do.call(cbind, equalized_dfs)
  
  return(combined_df)
}

#' Pad a Tibble to a Desired Length
#'
#' This function expands a tibble to a desired length by adding additional rows.
#' These new rows are filled with a specified value or `NA` by default.
#'
#' @param tbl A tibble or dataframe that you want to pad with additional rows.
#' @param desired_length The total number of rows you want the tibble to have after padding.
#' @param fill (optional) The value used to fill the new rows. By default, this is `NA`.
#'   The fill value must be compatible with the types of the columns in the tibble.
#'
#' @return A tibble that has been padded to the desired length.
#'
#' @examples
#' example_tbl <- tibble(Assay = c("ADA2", "ADAM8"), UniProt = c("Q9NZK5", "P78325"))
#' padded_tbl <- pad_tibble_to_length(example_tbl, desired_length = 5)
#' print(padded_tbl)
#'
#' @export
pad_tibble_to_length <- function(tbl, desired_length, fill = NA) {
  current_length <- nrow(tbl)
  additional_rows <- desired_length - current_length
  
  # Check if additional rows are needed
  if (additional_rows > 0) {
    # Create a tibble with the additional rows filled with NA
    pad_tbl <- replicate(ncol(tbl), rep(fill, additional_rows), simplify = FALSE) %>%
      as_tibble() %>%
      setNames(names(tbl))
    
    # Bind the original tibble with the padding tibble
    tbl <- bind_rows(tbl, pad_tbl)
  }
  
  tbl
}


#' Helper Functions for Manipulation of Bone Marrow Olink Data
#'
#' These functions help in common operations such as filtering
#'
#' @name data.bm.helpers
#' @export
#' @importFrom dplyr select filter
#' @importFrom magrittr %>%
data.bm.helpers <- list(
select_data_columns = function(df) {
  df %>% dplyr::select(Assay, AssayUnique, OlinkID, UniProt, Panel, Visit, Subject, NPX, SampleID, SampleType)
},

filter_non_Ig_depleted = function(df) {
  df %>% dplyr::filter(SampleType == "IG+") %>% dplyr::mutate(SampleType = forcats::fct_drop(SampleType))
},

filter_Ig_depleted = function(df) {
  df %>% dplyr::filter(SampleType == "IG-") %>% dplyr::mutate(SampleType = forcats::fct_drop(SampleType))
},

filter_MM_samples = function(df) {
  df %>% dplyr::filter(is_MM_sample == TRUE & ! (Visit ==  "(unknown)")) %>% 
    dplyr::mutate(
      Visit = forcats::fct_drop(Visit),
      SampleType = forcats::fct_drop(SampleType),
      Subject = forcats::fct_drop(Subject)
    )
},

sample_assays = function(df, n) {
  df %>% 
  dplyr::filter(AssayUnique %in% sample(AssayUnique, n)) 
}

)

#' Helper Functions for Manipulation of Bone Marrow Olink Data
#'
#' These functions help in common operations such as filtering
#'
#' @name data.plasma.helpers
#' @export
#' @importFrom dplyr select filter
#' @importFrom magrittr %>%
data.plasma.helpers <- 
  list(
    drop_unused_levels = function(df) {
      df %>% dplyr::mutate(across(where(is.factor), forcats::fct_drop))
    },
    sort_Subject = function(df) {
      df %>% dplyr::mutate(Subject = factor(Subject, levels = sort(levels(Subject))))
    },
    filter_FH1_samples = function(df) {
      df %>% 
        dplyr::filter(Cohort == "FH1") %>% 
        mutate(across(where(is.factor), forcats::fct_drop))
    },
    filter_MM_samples = function(df) {
      df %>% 
      dplyr::filter(Cohort == "FH1") %>% 
      dplyr::mutate(across(where(is.factor), forcats::fct_drop)) %>%
      dplyr::filter(visitDetails %in% c(
                                 "MM Dx",
                                 "MM C2",
                                 "MM C4",
                                 "Tx 60",
                                 "Tx 90",
                                 "Tx 1Y"
                                )
                 ) %>%
      dplyr::mutate(across(where(is.factor), forcats::fct_drop))
    }
)


