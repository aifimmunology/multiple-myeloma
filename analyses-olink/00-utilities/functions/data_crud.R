#' Load a single file
#'
#' This function loads a single CSV file based on a provided file handle. It reads the file into R as a data.frame or tibble,
#' ensuring that string data is not converted into factors. Errors during file reading are logged and reported.
#'
#' @param handle A list containing the file path under the $file key.
#' @importFrom logger log_info log_error
#' @importFrom readr read_csv
#' @return A data.frame or tibble from the read CSV file.
#' @export
#' @examples
#' handle <- list(file = "path/to/your/file.csv")
#' df <- loadFile(handle)
loadFile <- function(handle) {
  file_path <- handle$file
  if (!file.exists(file_path)) {
    stop(sprintf("Failed to read: File does not exist at path: '%s'", file_path))
  }
  tryCatch({
    logger::log_info(sprintf("Attempting to read file at '%s'", file_path))
    readr::read_csv(file_path, stringsAsFactors = FALSE)
  }, error = function(e) {
    error_msg <- sprintf("Error reading file at '%s': %s", file_path, e$message)
    logger::log_error(error_msg)
    stop(error_msg)
  })
}

#' Load multiple files
#'
#' This function reads multiple files specified by a list of data handles. Each data handle should contain a file path.
#' The function can optionally read a single file specified by its index within the list. It leverages `loadFile` for each file,
#' ensuring consistency in file reading behavior.
#'
#' @param data.handles A list of lists, each containing a file path under the $file key.
#' @param idx Optional index to specify a particular file to load.
#' @importFrom logger log_info
#' @export
#' @examples
#' data_handles <- list(list(file = "path/to/your/first_file.csv"), list(file = "path/to/your/second_file.csv"))
#' all_files <- loadFiles(data_handles)
#' specific_file <- loadFiles(data_handles, idx = 1)
loadFiles <- function(data.handles, idx = NULL) {
  stopifnot(is.list(data.handles), all(sapply(data.handles, function(x) is.list(x) && !is.null(x$file))))
  
  logger::log_info("Starting to read cache files")
  
  if (!is.null(idx)) {
    stopifnot(idx %in% seq_along(data.handles))
    file_contents <- list(loadFile(data.handles[[idx]]))
  } else {
    file_contents <- purrr::map(data.handles, loadFile)
  }
  
  return(file_contents)
}

#' Retrieve file names from HISE handles
#'
#' Extracts and returns the file names from a list of HISE data handles. This is useful for operations that require only
#' the file names from a set of HISE data handles, such as logging or further file manipulations.
#'
#' @param hise.handles A list of HISE data handles.
#' @importFrom purrr map_chr
#' @return A character vector of file names.
#' @export
#' @examples
#' hise_handles <- list(list(file = "first_file.csv"), list(file = "second_file.csv"))
#' file_names <- getFileNames(hise_handles)
getFileNames <- function(hise.handles) {
  purrr::map_chr(hise.handles, ~ .x$file)
}

#' Download files from HISE
#'
#' Downloads files from the HISE database using a list of file IDs. It logs the download process and returns a list of data handles
#' for the downloaded files. Each handle includes the local file path, facilitating subsequent file operations.
#'
#' @param ids A list of file IDs to download.
#' @importFrom logger log_info log_debug
#' @importFrom hise readFiles
#' @return A list of data handles for the downloaded files.
#' @export
#' @examples
#' file_ids <- list("file_id_1", "file_id_2")
#' downloaded_files <- downloadFiles(file_ids)
downloadFiles <- function(ids) {
  if (length(ids) == 0) {
    stop("The 'ids' argument must be a non-empty list of file IDs.")
  }

  logger::log_info("Downloading files from HISE")
  logger::log_debug("Attempting to retrieve the following ids: {ids}")
  data.handles.raw <- hise::readFiles(ids)
  
  if (is.null(data.handles.raw) || length(data.handles.raw) == 0) {
    logger::log_info("No files were downloaded from HISE.")
    return(NULL)
  } else {
    file.paths <- getFileNames(data.handles.raw)
    logger::log_info("Files downloaded from HISE and cached")
    logger::log_debug("Downloaded files: {file.paths}")
    return(data.handles.raw)
  }
}

#' Read CSV files into R
#'
#' Reads multiple CSV files specified by file paths into R as tibbles. It performs checks to ensure each file exists before
#' attempting to read, and logs each step of the process.
#'
#' @param file_names A character vector of file paths to read.
#' @importFrom readr read_csv
#' @importFrom logger log_info
#' @return A list of tibbles from the read CSV files.
#' @export
#' @examples
#' file_names <- c("path/to/your/first_file.csv", "path/to/your/second_file.csv")
#' csv_data <- readCSVFiles(file_names)
readCSVFiles <- function(file_names, ...) {
  if (!is.character(file_names) || length(file_names) == 0) {
    stop("file_names must be a non-empty character vector.", call. = FALSE)
  }
  
  missing_files <- file_names[!file.exists(file_names)]
  if (length(missing_files) > 0) {
    stop("One or more files do not exist.", call. = FALSE)
  }
  
  logger::log_info("Starting to read CSV files.")
  tibbles_list <- Map(function(file_name) {      #Map keeps names
    logger::log_info(paste("Reading file:", file_name))
    data_csv <- readr::read_csv(file_name, ...)
    
  },
  file_names)
  if(length(tibbles_list) == 0) warning("No files read. Returning empty list")
  logger::log_info("Completed reading all CSV files.")
  
  return(tibbles_list)
}

#' @describeIn readCSVFiles Read a single CSV file into R as a tibble.
#'
#' @inherit readCSVFiles
#' @export 
#' @examples
#' file_name <- "path/to/your/file.csv"
#' csv_data <- readCSVFile(file_name)
readCSVFile <- function(file_name, ...) {
  dfs.list <- readCSVFiles(file_name, ...)
  
  if(length(dfs.list) == 1) return(dfs.list[[1]])
  warning("More than on file read. Returning list")
  return(dfs.list)
}

                                             
#' From Lauren Okada. Extract Specimen-to-Kit ID Key from File Descriptors
#'
#' @param fd List type. A result returned by hise::getFileDescriptors()
#' @param keep_ids Character vector or list. Specimen ids you wish to keep, for
#' example the specimens in an Olink file "SampleID" column.
#' @param keep_fields Character vector, default c('specimenGuid','specimenType').
#' The fields from the 'specimens' file descriptors that you wish to keep.
#' @return A data.frame with rows as unique specimens, filtered to only keep_ids, if provided. 
#' Columns are 'sample.sampleKitGuid' and keep_fields. If no keep_ids, all specimens returned. If
#' no keep_fields, all fields returned.
#' @export
get_specimen_key <- function(fd, keep_ids = c(), keep_fields = c('specimenGuid','specimenType')){
    fd <- unlist(fd, recursive = FALSE) # fd are listed if multiple filter values used (only tested for 2 values of file id)

    dlist <- lapply(fd, function(x){
        spec <- lapply(x[['specimens']], as.data.frame) %>%
            bind_rows()
        spec$sample.sampleKitGuid <- unique(x$sample$sampleKitGuid) 
        spec
    })
    df <- dplyr::bind_rows(dlist)

    if(any(is.na(keep_ids)) || length(keep_ids)>0){
        keep_ids <- as.character(keep_ids)
        diff_id <- setdiff(keep_ids, df$specimenGuid)
        if(length(diff_id)>0){
            message(sprintf("%s specimen ids requested not found in file descriptors: %s",
                    length(diff_id),
                    paste(diff_id, collapse = ", ")))
        }
        match_ids <- intersect(keep_ids, df$specimenGuid)
        if(length(match_ids) >1){
            df <- df %>%
                dplyr::filter(specimenGuid %in% match_ids)
        } else {
            message(sprintf("No requested specimen id's found in file descriptors. Returning all specimens"))
        }
    } 

    if(any(is.na(keep_fields)) || length(keep_fields)>0){
        keep_fields <- unique(c(keep_fields, 'sample.sampleKitGuid'))
        diff_field <- setdiff(keep_fields, names(df))
        if(length(diff_field)>0){
            message(sprintf("%s fields requested not found in file descriptors: %s",
                    length(diff_field),
                    paste(diff_field, collapse = ", ")))
        }
        match_fields <- intersect(keep_fields, names(df))
        if(length(match_fields) >1){
            df <- df %>%
                dplyr::select(all_of(match_fields))
        } else {
            message(sprintf("None of the requested fields found in file descriptors. Returning all specimen fields"))
        }
    }

    df

}

#' Check if a File is in a Private Folder
#'
#' This function determines if a given file is located within a specified private folder by using a private API from the `hise` package. It checks all files in all private folders and verifies if the specified file exists in the specified folder.
#'
#' @param filename A character string representing the name of the file to check.
#' @param folder A character string representing the name of the folder to check the file against.
#'
#' @return A logical value (`TRUE` or `FALSE`) indicating whether the file is present in the specified private folder.
#'
#' @details The function uses the `listFilesInAllPrivateFolders` function from the `hise` package to retrieve a list of all files in all private folders. It then checks if the specified folder is among those folders and subsequently if the file exists in that folder. The function will stop and throw an error if the `filename` is not a character or if the specified `folder` does not exist among the private folders.
#'
#' @examples
#' # Assuming you have a hise package and appropriate private folders setup
#' # isInPrivateFolder("example.docx", "private_docs")
#'
#' @export
#' @importFrom hise listFilesInPrivateFolder
isInPrivateFolder <- function(folder, filename) {
  stopifnot(is.character(filename))
  stopifnot(is.character(folder))
  stopifnot(length(filename) == 1 & length(folder) == 1) #not implemented vectorized check yet
  hise.files.df <- hise::listFilesInPrivateFolder(folderName = folder, toDF = TRUE) %>% as_tibble
  
  folder.files <- hise.files.df %>%
    filter(folderName == folder) %>%
    pull(fileName)

  any(filename %in% folder.files)
  
}

                                              
#' Write a File to a Private Folder
#'
#' This function uploads a file to a specified private folder. If a file with the same name already exists in the folder, it is deleted before the new file is uploaded. This ensures that the file in the folder is always the most recent version. The function is designed to work with private folders specified by their unique identifiers.
#'
#' @param folder A character string representing the unique identifier of the private folder where the file will be uploaded.
#' @param filename A character string representing the name of the file to upload.
#'
#' @return Invisible NULL. The function is called for its side effects.
#'
#' @details The function first checks if the specified folder and filename are character strings and whether the folder is one of the allowed private folders. It uses the `isInPrivateFolder` function to check if a file with the same name already exists in the specified folder. If such a file exists, it logs the deletion of the old file using `logger::log_info` and then deletes the file with `hise::deleteFileInPrivateFolder`. Finally, it uploads the new file using `hise::uploadFileToPrivateFolder`.
#'
#' The allowed folders are currently hard-coded as "lt-ndmm" and "lt-ndmm-data".
#'
#' @examples
#' # Assuming you have configured the hise package, the logger package,
#' # and the specific private folders
#' # writePrivateFolder("lt-ndmm", "report.pdf")
#'
#' @export
#' @importFrom hise deleteFileInPrivateFolder uploadFileToPrivateFolder
#' @importFrom logger log_info
writePrivateFolder <- function(folder, file_path) {
  stopifnot(is.character(folder))
  stopifnot(is.character(file_path))
  stopifnot(folder %in% c("lt-ndmm", "lt-ndmm-data"))
  stopifnot(file.exists(file_path))
  if(isInPrivateFolder(folder, file_path)) {
    logger::log_info("Deleting old file of same name: {filename}")
    hise::deleteFileInPrivateFolder(folder, fileName = basename(file_path))
  }

  hise::uploadFileToPrivateFolder(folderName = folder, fileName = file_path)
  
}

checkDataFrame <- function(df) {
  # Check if the input is a dataframe
  if (!is.data.frame(df)) {
    message("The input is not a dataframe.")
    return(FALSE)
  }
  
  # Check if the dataframe has any data
  if (nrow(df) <= 0 || ncol(df) <= 0) {
    message("The dataframe is empty (no rows or no columns).")
    return(FALSE)
  } 

  TRUE
}


#' Write DataFrame to CSV in a Private Folder
#'
#' This function writes a given dataframe to a CSV file intended for Excel compatibility,
#' and then uploads the file to a specified private folder. It first checks that the inputs
#' for folder and filename are character strings and that the folder is one of the allowed
#' private folders. It also verifies that the dataframe is a valid dataframe and contains data
#' before proceeding to write the file. The CSV file is then uploaded to the private folder.
#'
#' @param df A dataframe to be written to a CSV file. The dataframe should not be empty
#' and must be a valid R dataframe object.
#' @param folder A character string specifying the private folder where the CSV file will be uploaded.
#' The folder must be one of the predefined allowed folders ("lt-ndmm" or "lt-ndmm-data").
#' @param filename A character string specifying the name of the CSV file to be created and uploaded.
#' The filename should include the `.csv` extension.
#' @param ... Additional arguments to be passed to `readr::write_excel_csv`, allowing customization
#' of the CSV file format and write options.
#'
#' @return Invisible NULL. The function is called for its side effects of writing a file to disk
#' and uploading it to a private folder.
#'
#' @details The function uses `checkDataFrame` to verify the dataframe's validity and contents.
#' It employs `readr::write_excel_csv` for creating a CSV file that is compatible with Excel.
#' After writing the file to disk, it uses `writePrivateFolder` to upload the file to the specified
#' private folder. Note that this function assumes `checkDataFrame` and `writePrivateFolder`
#' are defined and available in the environment.
#'
#' @examples
#' \dontrun{
#' df <- data.frame(x = 1:5, y = letters[1:5])
#' writeCSVPrivateFolder(df, "lt-ndmm", "example.csv")
#' }
#'
#' @importFrom readr write_excel_csv
#' @export
writeCSVPrivateFolder <- function(df, folder, file_path, ...) {
  stopifnot(is.character(folder))
  stopifnot(is.character(file_path))
  stopifnot(folder %in% c("lt-ndmm", "lt-ndmm-data"))
  
  # Ensure df is a valid dataframe with data
  stopifnot(checkDataFrame(df))

  ensure_folder(file_path)
  # Write the dataframe to a CSV file for Excel compatibility
  readr::write_excel_csv(df, file_path, ...)
  
  if(isInPrivateFolder(folder, basename(file_path))) {
    hise::deleteFileInPrivateFolder(folder, basename(file_path))
  }
  
  writePrivateFolder(folder, file_path)
}

ensure_folder <- function(file_path) {
    directory = dirname(file_path)
    if(! dir.exists(directory)) {
        dir.create(directory, recursive = TRUE)
    }
}
                                              
#' @export
writeRDSPrivateFolder <- function(object, folder, file_path) {
  ensure_folder(file_path)
    
  saveRDS(object, file_path)
  
  if(isInPrivateFolder(folder, basename(file_path))) {
      hise::deleteFileInPrivateFolder(folder, basename(file_path))
  }

  writePrivateFolder(folder, file_path)
  
}
                                              
#' Workhorse for Downloading From Private Folder
#' 
#' Downloads files and reads them with read_csv
#'
#' @return a tibble
#' @importFrom readr read_csv
downloadPrivateFolder <- function(folder, filename, ...) {
  if(! isInPrivateFolder(folder, filename)) {
    stop(paste("The file", filename, "could not be found in", folder))
  }
  hise::downloadFileFromPrivateFolder(folderName = folder, fileName = filename)
  
  df <- readr::read_csv(filename, ...)
  if(file.exists(filename)) {
      file.remove(filename)
  }
  df
}

                                              
#' Retrieve Cleaned FH1 Assembled Labs Dataset
#'
#' Downloads the FH1 Assembled Labs dataset from a specified private folder and reads it into R.
#' The dataset is expected to be in CSV format, and a predefined column specification is used to ensure
#' correct data types are applied upon import.
#'
#' @return A dataframe containing the FH1 Assembled Labs dataset.
#' @importFrom hise downloadFileFromPrivateFolder
#' @importFrom readr read_csv
#' @export
#' @examples
#' \dontrun{
#' df <- getFH1AssembledLabs()
#' }
getFH1AssembledLabs <- function() {
  filename = "FH1 Assembled Labs w April Data - cleaned.csv"
  folder = "lt-ndmm"
  col_types = FH1AssembledLabsColumnSpec
  downloadPrivateFolder(folder= folder, filename = filename, col_types = col_types)
}

#' Retrieve Cleaned FH1 Treatment Data
#'
#' Downloads Treatment Data for FH1
#'
#' @return A dataframe 
#' @importFrom readr read_csv
#' @export
getFH1Treatment <- function() {
  filename = "FH1 Treatment Cycles Feb 2024 - cleaned.csv"
  folder = "lt-ndmm"
  downloadPrivateFolder(folder = folder, filename = filename)
}

#' Retrieve MM vs Healthy Stat 
#'
#' Downloads Stats from Comparison of Healthy BR2 vs MM
#'
#' @return A dataframe 
#' @importFrom readr read_csv
#' @export
getLMERMPlasmaMvsHealthyBR2 <- function() {
  filename = "MM vs Healthy BR2 stats.csv"
  folder = "lt-ndmm"
  df <- downloadPrivateFolder(folder = folder, filename = filename)
  df %>%
    dplyr::mutate(contrast =  forcats::fct(contrast, 
                                           levels = c('MM Dx - Healthy', 'MM C2 - Healthy', 'MM C4 - Healthy', 'Tx 90 - Healthy', 'Tx 1Y - Healthy')
                                          ),
                 MM.status = forcats::fct(MM.status, 
                                          levels = c('MM Dx', 'MM C2', 'MM C4', 'Tx 90', 'Tx 1Y')
                                         )
                 )
}

#' Retrieve MM vs Healthy Stat
#'
#' Downloads Stats from Comparison of Healthy BR1 and BR2 vs MM
#'
#' @return A dataframe 
#' @importFrom readr read_csv
#' @export
getLMERPlasmaMMvsHealthy <- function() {
  filename = "MM vs Healthy BR1 and BR2 stats.csv"
  folder = "lt-ndmm"
  df <- downloadPrivateFolder(folder = folder, filename = filename)
  df %>%
    dplyr::mutate(contrast =  forcats::fct(contrast, 
                                           levels = c('MM Dx - Healthy', 'MM C2 - Healthy', 'MM C4 - Healthy', 'Tx 90 - Healthy', 'Tx 1Y - Healthy')
                                          ),
                 MM.status = forcats::fct(MM.status, 
                                          levels = c('MM Dx', 'MM C2', 'MM C4', 'Tx 90', 'Tx 1Y')
                                         )
                 )
}
                                              
#' Retrieve Cleaned Bone Marrow Olink Data
#'
#' Downloads Olink Data for FH1
#'
#' @return A dataframe 
#' @importFrom readr read_csv
#' @export
getDataOlinkBM <- function() {
  filename = "Data Olink BM - cleaned.csv"
  folder = "lt-ndmm"
  col_types = DataBMColumnSpec
  df <- downloadPrivateFolder(folder, filename, col_types = col_types)

  df %>%
    dplyr::mutate(Visit = 
           forcats::fct_relevel(Visit, c("MM Dx", "MM C4", "MM Tx90", "MM Tx1Y", "MM Tx2Y", "Normal BM", "(unknown)"))
           )
           
}

#' Retrieve Cleaned Plasma Olink Data
#'
#' Downloads Olink Data for FH1, BR1, BR2 Plasma
#'
#' @return A dataframe 
#' @importFrom readr read_csv
#' @export
getDataOlinkPlasma <- function() {
  filename = "Olink_Plasma_FH1_BR1_BR2_with_hise_descriptors.csv"
  folder = "lt-ndmm"
  col_types = DataPlasmaColumnSpec
  df <- downloadPrivateFolder(folder, filename, col_types = col_types)
  df %>%
  dplyr::mutate(
    visitDetails = forcats::fct_recode(visitDetails,
                     `MM Dx` = "MM Pre-Treatment",
                     `MM C2` = "MM Post Induction 2-Cycles",
                     `MM C4` = "MM End Induction 1st Draw",
                     `Tx 60` = "MM Post Transplant 60 Days",
                     `Tx 90` = "MM Post Transplant 90 Days",
                     `Tx 1Y` = "MM Post Transplant 1 year",                     
                     `Flu` = "N/A - Flu-Series Timepoint Only",
                     `Standalone` = "N/A - stand-alone collection",
                     `Other` = "Other",
                     `Imm d0` = "Immune Variation Day 0",
                     `Imm d7` = "Immune Variation Day 7",
                     `Imm d90` = "Immune Variation Day 90"
                ), #ensure order
      visitDetails = forcats::fct_relevel(visitDetails,
                                 "MM Dx",
                                 "MM C2",
                                 "MM C4",
                                 "Tx 60",
                                 "Tx 90",
                                 "Tx 1Y",
                                 "Flu",
                                 "Standalone",
                                 "Other",
                                 "Imm d0",
                                 "Imm d7",
                                 "Imm d90"
                )
        ) %>%
    dplyr::mutate(Age = lubridate::year(drawDate) - birthYear) %>%
    dplyr::mutate( #reorder Subjects according to two last digits in subject identifier
        Subject = forcats::fct_reorder(Subject, 
                              as.numeric(
                                  substring(
                                      as.character(Subject), 
                                      nchar(as.character(Subject))-1, 
                                      nchar(as.character(Subject))
                                  )
                              )
                             )
    )
}

#' Retrieve Stats Results from LMER Mixed Model Analysis from Bone Marrow
#'
#' Downloads stats from lmer analysis on bone amrrow samples
#'
#' @return A dataframe 
#' @importFrom readr read_csv
#' @export                                              
getLMERSignificantProteinsBM <- function() {
  filename = "Stats Proteins lmer by visit Bone Marrow Non IG Depleted.csv"
  folder = "lt-ndmm"
  df.stats <- downloadPrivateFolder(folder, filename)


}

#' Retrieve Stats Results from LMER Mixed Model Analysis from Plasma
#'
#' Downloads stats from lmer analysis on plasma samples
#'
#' @return A dataframe 
#' @importFrom readr read_csv
#' @export                                                
getLMERSignificantProteinsPlasma <- function() {
  filename = "Stats Proteins lmer by visit Plasma.csv"
  folder = "lt-ndmm"
  df.stats <- downloadPrivateFolder(folder, filename)
#  data.bm <- getDataOlinkBM()
  df.stats %>%
      dplyr::mutate(contrast = fct(contrast, levels = c("MM C2 - MM Dx", "MM C4 - MM Dx", "MM C4 - MM C2", "Tx 60 - MM C4", "Tx 1Y - Tx 60")))
}


#' Write Data Frame to JSON in a Private Folder
#'
#' This function checks if the directory for the provided file path exists, creates it if it doesn't,
#' then writes the given data frame to a JSON file at the specified path. If the file is in a private
#' folder, it is first deleted using hise::deleteFileInPrivateFolder before being written with 
#' ndmmFH1::writePrivateFolder.
#'
#' @param df Data frame to be written as JSON.
#' @param folder The name or path of the folder considered as private.
#' @param file_path The full path where the JSON file will be saved. This includes the directory
#'   path and the file name.
#' @importFrom jsonlite write_json
#' @importFrom hise deleteFileInPrivateFolder
#' @examples
#' \dontrun{
#'   df <- data.frame(name = c("Alice", "Bob"), age = c(30, 32))
#'   writeJSON(df, "myPrivateFolder", "path/to/myPrivateFolder/data.json")
#' }
#' @export
writeJSON <- function(df, folder, file_path) {
    dir.path <- dirname(file_path)
    if (!dir.exists(dir.path)) {
      dir.create(dir.path, recursive = TRUE)
    }
    
    #js <- jsonlite::toJSON(df, dataframe = "rows", pretty = TRUE)
    jsonlite::write_json(df, dataframe = "rows", pretty = TRUE, path = file_path)
    
    if(isInPrivateFolder(folder, basename(file_path))) {
      hise::deleteFileInPrivateFolder(folder, basename(file_path))
    }

    writePrivateFolder(folder, file_path)
}
