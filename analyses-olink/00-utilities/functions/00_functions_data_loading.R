library(hise)
require(data.table)
require(ggplot2)
require(dplyr)
library(patchwork)
library(ggpubr)
library(tidyverse)
library(magrittr)

# Love Tatting ndmmFH1/R/static.R
DataBMColumnSpec <- readr::cols(
  SampleID = readr::col_factor(),
  Index = readr::col_factor(),
  OlinkID = readr::col_character(),
  UniProt = readr::col_character(),
  Assay = readr::col_character(),
  MissingFreq = readr::col_double(),
  Panel = readr::col_factor(),
  Panel_Lot_Nr = readr::col_factor(),
  PlateID = readr::col_factor(),
  QC_Warning = readr::col_factor(),
  LOD = readr::col_double(),
  NPX = readr::col_double(),
  Normalization = readr::col_factor(),
  Assay_Warning = readr::col_factor(),
  is_technical_control = readr::col_logical(),
  SampleType = readr::col_factor(),
  Subject = readr::col_factor(),
  Sex = readr::col_factor(),
  Age = readr::col_double(),
  Race = readr::col_factor(),
  Cohort = readr::col_factor(),
  Visit = readr::col_factor(),
  CollectionDate = readr::col_date(format = "%Y-%m-%d"),
  sampleKitGuid = readr::col_character(),
  is_normal_BM = readr::col_logical(),
  is_MM_sample = readr::col_logical(),
  AssayUnique = readr::col_factor()
)

DataPlasmaColumnSpec <- readr::cols(
  SampleID = readr::col_character(),
  Index = readr::col_double(),
  AssayUnique = readr::col_character(),
  OlinkID = readr::col_character(),
  UniProt = readr::col_character(),
  Assay = readr::col_character(),
  MissingFreq = readr::col_double(),
  Panel = readr::col_character(),
  Panel_Lot_Nr = readr::col_factor(),
  PlateID = readr::col_factor(),
  Normalization = readr::col_character(),
  Assay_Warning = readr::col_character(),
  sampleKitGuid = readr::col_character(),
  UniqueAssay = readr::col_character(),
  N_Bridge = readr::col_factor(),
  BatchOffset = readr::col_double(),
  NPX = readr::col_double(),
  LOD = readr::col_double(),
  sampleWarningBridge = readr::col_character(),
  BridgingDetails = readr::col_character(),
  filename = readr::col_character(),
  BatchID = readr::col_factor(),
  is_bridgingControl = readr::col_double(),
  visitName = readr::col_factor(),
  visitDetails = readr::col_factor(),
  drawDate = readr::col_datetime(format = "%Y-%m-%d"),
  daysSinceFirstVisit = readr::col_double(),
  diseaseStatesRecordedAtVisit = readr::col_factor(),
  file.fileType = readr::col_character(),
  Sex = readr::col_factor(),
  birthYear = readr::col_double(),
  Ethnicity = readr::col_factor(),
  Race = readr::col_factor(),
  Subject = readr::col_factor(),
  Cohort = readr::col_factor(),
  sampleWarningBride = readr::col_factor()
)

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
  hise.files.df <- hise::listFilesInPrivateFolder(
    folderName = folder,
    toDF = TRUE
  ) %>%
    as_tibble

  folder.files <- hise.files.df %>%
    filter(folderName == folder) %>%
    pull(fileName)

  any(filename %in% folder.files)
}

downloadPrivateFolder <- function(folder, filename, ...) {
  if (!isInPrivateFolder(folder, filename)) {
    stop(paste("The file", filename, "could not be found in", folder))
  }
  hise::downloadFileFromPrivateFolder(folderName = folder, fileName = filename)

  df <- readr::read_csv(filename, ...)
  if (file.exists(filename)) {
    file.remove(filename)
  }
  df
}

filter_FH1_samples = function(df) {
  df %>%
    dplyr::filter(Cohort == "FH1") %>%
    mutate(across(where(is.factor), forcats::fct_drop))
}

# Love Tatting ndmmFH1::getDataOlinkPlasma
getDataOlinkPlasma <- function() {
  filename = "Olink_Plasma_FH1_BR1_BR2_with_hise_descriptors.csv"
  folder = "lt-ndmm"
  col_types = DataPlasmaColumnSpec
  df <- downloadPrivateFolder(folder, filename, col_types = col_types)
  df %>%
    dplyr::mutate(
      visitDetails = forcats::fct_recode(
        visitDetails,
        `MM Dx` = "MM Pre-Treatment",
        `MM C2` = "MM Post Induction 2-Cycles",
        `MM C4` = "MM End Induction 1st Draw",
        `Tx 60` = "MM Post Transplant 60 Days",
        `Tx 90` = "MM Post Transplant 90 Days",
        `Tx 1Y` = "MM Post Transplant 1 year",
        Flu = "N/A - Flu-Series Timepoint Only", # flu day 0, 7, 90
        Standalone = "N/A - stand-alone collection", # sometimes flu day 0, sometimes before day 0
        Other = "Other",
        `Imm d0` = "Immune Variation Day 0",
        `Imm d7` = "Immune Variation Day 7",
        `Imm d90` = "Immune Variation Day 90"
      ),
      visitDetails = forcats::fct_relevel(
        visitDetails,
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
    dplyr::mutate(
      Subject = forcats::fct_reorder(
        Subject,
        as.numeric(substring(
          as.character(Subject),
          nchar(as.character(Subject)) -
            1,
          nchar(as.character(Subject))
        ))
      )
    )
}

# Love Tatting ndmmFH1::getDataOlinkBM
getDataOlinkBM <- function() {
  filename = "Data Olink BM - cleaned.csv"
  folder = "lt-ndmm"
  col_types = DataBMColumnSpec
  df <- downloadPrivateFolder(folder, filename, col_types = col_types)
  df %>%
    dplyr::mutate(
      Visit = forcats::fct_relevel(
        Visit,
        c(
          "MM Dx",
          "MM C4",
          "MM Tx90",
          "MM Tx1Y",
          "MM Tx2Y",
          "Normal BM",
          "(unknown)"
        )
      )
    )
}

annotateFluResponse <- function(df) {
  # Define the flu_list
  Responders <- c(
    'FH1002',
    'FH1005',
    'FH1006',
    'FH1008',
    'FH1012',
    'FH1014',
    'FH1017'
  )
  NonResponders <- c(
    'FH1021',
    'FH1016',
    'FH1011',
    'FH1009',
    'FH1007',
    'FH1004',
    'FH1003'
  )

  df$manual.flu_response <- "N/A"
  df <- df %>%
    dplyr::mutate(
      manual.flu_response = case_when(
        Subject %in% Responders ~ "responder",
        Subject %in% NonResponders ~ "non_responder"
      )
    )

  df
}

load_data_plasma_flu <- function(flu_series = FALSE) {

  fluVisitNames <- c(
      'Flu Year 1 Day 0',
      'Flu Year 1 Day 7',
      'Flu Year 1 Day 90',
      'Flu Year 2 Day 0',
      'Flu Year 2 Day 7',
      'Flu Year 2 Day 90'
  )

  plasma <- setDT(getDataOlinkPlasma())

  # Rename standalone to Day 0
  # Justification: standalone collection is all 
  # on day 0 or prior
  plasma[visitName == 'Flu Year 1 Stand-Alone', visitName := "Flu Year 1 Day 0"]
  plasma[visitName == 'Flu Year 2 Stand-Alone', visitName := "Flu Year 2 Day 0"]

  # Filter to these visitNames
  # Here we do not filter to just FH1 cohort
  # BR1 and BR2 cohorts have the same fluVisitNames :))
  plasma <- plasma %>%
    filter(visitName %in% fluVisitNames)
  plasma <- setDT(plasma)

  # ONLY analyse VRD subjects, remove FH1 subjects treated with dara
  dara_pos <- c(
    'FH1020',
    'FH1022',
    'FH1023',
    'FH1024',
    'FH1026',
    'FH1027',
    'FH1028'
  )
  plasma[, dara := Subject %in% dara_pos]
  plasma = plasma[dara == FALSE]

  # Annotate flu response
  plasma <- annotateFluResponse(plasma)

  # For downstream compatibility, overwrite the 
  # visitDetails column to the ordered visitName for flu
  plasma$visitDetails_old <- plasma$visitDetails
  plasma$visitDetails <- factor(
    plasma$visitName,
    levels = fluVisitNames
  )

  plasma <- plasma[, c(
    'NPX',
    'Assay',
    'visitDetails',
    'visitDetails_old',
    'visitName',
    'Subject',
    'SampleID',
    'Age',
    'Sex',
    'Cohort',
    'daysSinceFirstVisit',
    'drawDate',
    'dara'
  )]

  plasma
}

load_data_plasma <- function() {
  plasmaVisitNames <- c(
    "MM Dx",
    # "MM C2", commented out unused timepoints
    "MM C4",
    "Tx 60",
    # "Tx 90", # 2611 assays
    "Tx 1Y",
    "MM Post Transplant 2 year" # 2580 assays
  )

  plasma <- getDataOlinkPlasma()
  plasma <- plasma %>%
    filter_FH1_samples() %>%
    filter(visitDetails %in% plasmaVisitNames)
  plasma <- setDT(plasma)

  dara_pos <- c(
    'FH1020',
    'FH1022',
    'FH1023',
    'FH1024',
    'FH1026',
    'FH1027',
    'FH1028'
  )

  # ONLY analyse VRD subjects, remove subjects treated with dara
  plasma[, dara := Subject %in% dara_pos]
  plasma = plasma[dara == FALSE]
  plasma <- plasma[, c(
    'NPX',
    'visitDetails',
    'visitName',
    'Subject',
    'SampleID',
    'Assay',
    'dara'
  )]

  # Annotate flu response
  plasma <- annotateFluResponse(plasma)

  # Rename Plasma 2Y timepoint
  plasma[visitDetails == 'MM Post Transplant 2 year', visitDetails := "Tx 2Y"]
  ### Set factors by visit
  plasma$visitDetails = factor(
    plasma$visitDetails,
    levels = c('MM Dx', 'MM C4', 'Tx 60', 'Tx 1Y', "Tx 2Y")
  )

  plasma
}

load_data_bm <- function(depleted = FALSE) {
  bmVisitNames <- c(
    "MM Dx",
    "MM C4",
    "MM Tx90",
    "MM Tx1Y",
    "MM Tx2Y"
  )

  bm <- getDataOlinkBM()

  # BM only, rename Visit column
  bm$visitDetails <- bm$Visit

  if (depleted == TRUE) {
    # BM only, remove IG depleted
    bm <- bm %>% dplyr::filter(grepl("-001-001", bm$SampleID))
    all(bm$SampleType == "IG-") # Must be TRUE
  } else {
    bm <- bm %>% dplyr::filter(grepl("-001-003", bm$SampleID))
    all(bm$SampleType == "IG+") # Must be TRUE
  }

  bm <- bm %>% filter_FH1_samples() %>% filter(visitDetails %in% bmVisitNames)
  bm <- setDT(bm)

  dara_pos <- c(
    'FH1020',
    'FH1022',
    'FH1023',
    'FH1024',
    'FH1026',
    'FH1027',
    'FH1028'
  )

  # ONLY analyse VRD subjects, remove subjects treated with dara
  bm[, dara := Subject %in% dara_pos]
  bm = bm[dara == FALSE]

  bm <- bm[, c('NPX', 'visitDetails', 'Subject', 'SampleID', 'Assay', 'dara')]

  bm <- annotateFluResponse(bm)

  bm$visitDetails = factor(
    bm$visitDetails,
    levels = c('MM Dx', 'MM C4', 'MM Tx90', 'MM Tx1Y', 'MM Tx2Y')
  )

  bm
}

plot_longitudinal_protein <- function(assay, dat) {
  df = dat[Assay == assay]
  df[, NPX_combined := mean(NPX), by = list(Subject, visitDetails)]

  ggplot(
    df,
    aes(
      x = visitDetails,
      y = NPX_combined,
      group = Subject,
      col = Subject
    )
  ) +

    geom_boxplot(
      data = df,
      aes(x = visitDetails, y = NPX_combined, group = visitDetails),
      alpha = 0.4
    ) +
    geom_point() +
    geom_line(alpha = 0.4) +
    ggtitle(assay) +
    theme_minimal() +
    xlab('Visit') +
    ylab('NPX') +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1))
}

# https://github.com/aifimmunology/NDMM-Olink/blob/420708511fa4d19569327e26a7a86ca907811f87/ndmmFH1/R/data_crud.R#L70
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

#' from https://github.com/aifimmunology/olink-qc-report/blob/imran-newformat/R/format_olink_metadata.R
#' Format metadata for a given HISE fileid
#'
#' hise package version 2.14.0. Downloads file if not found in current directory
#' cache, reads the file to identify specimen ids.
#'
#' @param hise_fileid HISE File id of Olink file
#' @param keep_values A named vector of metadata to retain,
#' default c(SampleID = "specimen.specimenGuid", Cohort = "cohort.cohortGuid",
#' Subject = "subject.subjectGuid", Sex = "subject.biologicalSex",
#' VisitName = "sample.visitName", VisitDetails = "sample.visitDetails",
#' DaysSinceFirstVisit = "sample.daysSinceFirstVisit",
#' SampleKitGuid = "sample.sampleKitGuid")
#' @param sep Optional delimiter character of the Olink file, used if it is a CSV. Default is a comma.
#'
#' @returns A dataframe of sample metadata containing all samples in the olink
#' file and associated \code{keep_values} metadata columns.
#' @export
format_olink_metadata <- function(
    hise_fileid,
    keep_values = c(
      SampleID = "specimen.specimenGuid",
      Cohort = "cohort.cohortGuid",
      Subject = "subject.subjectGuid",
      Sex = "subject.biologicalSex",
      VisitName = "sample.visitName",
      VisitDetails = "sample.visitDetails",
      DaysSinceFirstVisit = "sample.daysSinceFirstVisit",
      SampleKitGuid = "sample.sampleKitGuid"
    ),
    alt_fileid = NULL,
    alt_fileType = NULL,
    fileType = "Olink",
    sep=","
) {
  # Location of files in new NextGen HISE IDEs
  file_list <- list.files(file.path("../../data/olink"), recursive = TRUE, full.names = TRUE, pattern = ".*[.]xlsx$|.*[.]csv$")
  fname <- file_list[grepl(hise_fileid, file_list)]
  
  if (length(fname) < 1) {
    invisible(hise::cacheFiles(fileIds = as.list(hise_fileid)))
  }
  file_list <- list.files(file.path("../../data/olink"), recursive = TRUE, full.names = TRUE, pattern = ".*[.]xlsx$|.*[.]csv$")
  fname <- file_list[grepl(hise_fileid, file_list)]

  if (tools::file_ext(fname) == "csv") {
    dat <- read.csv(fname, na = "NA", sep=sep)
  } else { 
    dat <- readxl::read_excel(fname, na = "NA")
  }
  all_specimens <- unique(dat$SampleID)
  cat(sprintf(
    "%s unique specimen IDs in file",
    length(all_specimens)
  ), sep = "\n")
  rm(dat)

  if (!is.null(alt_fileid) & !is.null(alt_fileType)) {
    fd <- hise::getFileDescriptors(fileType = alt_fileType, filter = list(file.id = alt_fileid))[[1]]
  } else {
    fd <- hise::getFileDescriptors(fileType = fileType, filter = list(file.id = hise_fileid))[[1]]
  }
  
  fd1 <- lapply(fd, function(x) {
    x[setdiff(names(x), c("lab", "survey", "emr", "surveyScheme", "specimens"))]
  })
  fd1 <- dplyr::bind_rows(lapply(fd1, function(x) {
    as.data.frame(t(unlist(x)))
  }))

  specimen_key <- format_specimen_list(all_specimens, fd)
  fd2 <- fd1 %>%
    dplyr::left_join(specimen_key, c("sample.sampleKitGuid"))

  # remove kits in fd that did not have an associated sample id in file
  i_extra_kits <- which(is.na(fd2$specimen.specimenGuid))
  if (length(i_extra_kits) > 0) {
    extra_kits <- paste(fd2$sample.sampleKitGuid[i_extra_kits], collapse = ",")
    fd2 <- fd2[-c(i_extra_kits), ]
    cat(sprintf(
      "Found %s kit IDs associated with HISE file metadata with no matching sample in file: %s",
      length(i_extra_kits),
      extra_kits
    ), sep = "\n")
  }

  fd2 <- fd2 %>%
    dplyr::select(all_of(keep_values))


  missing_id <- setdiff(all_specimens, fd2$SampleID)
  missing_string <- ifelse(length(missing_id) > 0,
    paste0(": [", paste(missing_id, collapse = ","), "]"),
    ""
  )
  cat(sprintf(
    "Found %s specimen IDs that match file\nNo metadata for %s specimen IDs in file%s",
    nrow(fd2),
    length(missing_id),
    missing_string
  ), sep = "\n")

  return(fd2)
}

# https://github.com/aifimmunology/olink-qc-report/blob/imran-newformat/R/format_specimen_list.R
#' Format HISE specimen list helper
#'
#' Formats a key of matching kit id's for specific specimens id's based on
#' \code{hise::getFileDescriptors()} result
#'
#' @param keep_specimen_ids Identified specimen id's to keep
#' @param fileDescriptors The getFileDescriptors() result call on the associated
#' file.
#' @returns A data.frame of 'sample.sampleKitGuid' and 'specimen.specimenGuid'
#' for each value of \code{keep_specimen_ids}
#' @export
format_specimen_list <- function(keep_specimen_ids, fileDescriptors) {
  specimen_list <- lapply(fileDescriptors, function(x) {
    x[["specimens"]]
  })
  kit_id <- sapply(fileDescriptors, function(x) {
    x[["sample"]]$sampleKitGuid
  })

  specimen_list <- lapply(specimen_list, function(x) {
    spec <- unlist(sapply(x, "[[", "specimenGuid"))
    spec[spec %in% keep_specimen_ids]
  })
  names(specimen_list) <- kit_id

  specimen_df_list <- lapply(kit_id, function(x) {
    if (is.null(specimen_list[[x]]) | length(specimen_list[[x]])==0) {
      data.frame(
        specimen.specimenGuid = NA,
        sample.sampleKitGuid = x
      )
    } else {
      data.frame(
        specimen.specimenGuid = specimen_list[[x]],
        sample.sampleKitGuid = x
      )
    }
  })
  specimen_key <- do.call(rbind, specimen_df_list)
  specimen_key <- specimen_key[!is.na(specimen_key$specimen.specimenGuid), ]

  return(specimen_key)
}