# ---
# jupyter:
#   jupytext:
#     formats: R:percent
#     text_representation:
#       extension: .R
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.3
#   kernelspec:
#     display_name: R (seurat_v5)
#     language: R
#     name: seurat_v5
# ---

# %%
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(stringr)
  library(tibble) 
})


# %%
#' Run Gene Set Enrichment Analysis (GSEA) per cell type
#'
#' This function performs Gene Set Enrichment Analysis (GSEA) on a list of genes for each cell type.
#'
#' @param logfc_list A data frame containing the gene expression data with columns for cell type, log fold change, gene, and other optional columns.
#' @param gmx A gene set matrix (GMT) file or a pre-loaded gene set matrix. If not provided, the function will load and initialize the pathway database using a default GMT file.
#' @param ct.col The column name in the logfc_list data frame that represents the cell type.
#' @param rank.col The column name in the logfc_list data frame that represents the log fold change.
#' @param gene.col The column name in the logfc_list data frame that represents the gene.
#' @param collapsePathways Logical value indicating whether to collapse pathways with similar gene sets into a single pathway.
#' @param ncores The number of cores to use for parallelization. If not provided, the function will use all available cores minus 3.
#'
#' @return A data frame containing the results of the GSEA analysis, including pathway information, p-values, adjusted p-values, normalized enrichment scores (NES), and leading edge genes.
#'
#' @import fgsea
#' @import dplyr
#' @import rlist
#' @import BiocParallel
#' @import data.table
#'
#' @examples
#' # Example usage of RunGSEACelltype function
#' logfc_list <- data.frame(cell_type = c("A", "A", "B", "B"),
#'                          logfc = c(1.2, 0.8, 1.5, 0.5),
#'                          gene = c("gene1", "gene2", "gene3", "gene4"))
#' result <- RunGSEACelltype(logfc_list)
#' print(result)
#'
#' @export

RunGSEACelltype <- function(logfc_list, gmx = NULL, ct.col = "cell_type", rank.col = "logfc", gene.col='gene',
                            collapsePathways = FALSE,
                            ncores = 1) {
    require(fgsea)
    # if no provided, Load and initialize pathway database
    if (is.null(gmx)) {
        gmxFile <- "../../../data/gmt/c2.cp.v7.2.symbols.gmt"
        colNames <- max(count.fields(file = gmxFile, sep = "\t"))
        colNames <- seq(from = 1, to = colNames)
        colNames <- as.character(colNames)
        gmx <- read.table(
            file = gmxFile,
            sep = "\t",
            quote = "\"",
            fill = TRUE,
            col.names = colNames,
            row.names = 1
        )
        gmx <- gmx[, -1]
        gmx <- apply(gmx, MARGIN = 1, FUN = function(x) {
            return(value = setdiff(unname(x), ""))
        })
        names(gmx) <- toupper(names(gmx))
    }

    # setup parallelization parameters
    if (is.null(ncores)) {
        ncores <- parallel::detectCores() - 3
    } else {
        (ncores <- ncores)
    }
    param <- BiocParallel::MulticoreParam(workers = ncores, progressbar = TRUE)

    # RUN GSEA per celltype
    celltypes <- unique(logfc_list %>% pull(.data[[ct.col]]))

    pLS <- lapply(celltypes, function(ct) {
        message(paste("run GSEA in", ct))

        # create rank list based on lowest to higest gene fold-change
        rnkDF <- logfc_list %>%
            dplyr::filter(.data[[ct.col]] == ct) %>%
            filter(!is.na(.data[[rank.col]])) %>%
            dplyr::arrange(.data[[rank.col]])
        rnk <- rnkDF %>%
            pull(.data[[rank.col]]) %>%
            as.numeric()
        names(rnk) <- rnkDF %>% pull(.data[[gene.col]])
        message(paste("run GSEA in", length(rnk), "genes"))

        # run GSEA by parallelization
        fgseaRes <- fgsea::fgsea(
            pathways = gmx,
            stats = rnk,
            minSize = 10,
            maxSize = 500, nPermSimple = 10000,
            BPPARAM = param
        )

        # filter on pathways <0.05 adjusted p-value
        fgseaRes_tb <- fgseaRes %>%
            as.data.frame() %>%
            # dplyr::filter(padj < 0.05) %>%
            dplyr::select(pathway, pval, padj, NES, leadingEdge) %>%
            dplyr::arrange(desc(NES)) %>%
            dplyr::mutate(celltype = ct)
        # if only keep the main pathway
        if (collapsePathways) {
            collapsedPathways <- fgsea::collapsePathways(
                fgseaRes[order(pval)][padj < 0.05],
                gmx, rnk
            )
            mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][
                order(-NES), pathway
            ]
            fgseaRes_tb <- fgseaRes_tb %>% dplyr::filter(pathway %in% mainPathways)
        }

        return(value = fgseaRes_tb)
    })

    pathwayDF <- data.table::rbindlist(pLS)
    pathwayDF$leadingEdge <- vapply(pathwayDF$leadingEdge,
        paste,
        collapse = ", ",
        character(1L)
    )

    # make plotting data frame
    plotDF <- pathwayDF %>%
        mutate(
            group = ifelse(NES > 0, "up", "down"),
            pID = c(1:length(pathway))
        ) %>%
        group_by(pID) %>%
        mutate(lesize = length(unlist(strsplit(leadingEdge, ",")))) %>%
        as_tibble()


    # determine pathway size
    gsSize <- data.frame(gsize = sapply(gmx, function(x) length(x))) %>%
        rownames_to_column(var = "pathway")

    # calculate propotion of genes enriched (#leading edge genes/size of pathway)
    plotDF <- plotDF %>%
        mutate(
            gsize = gsSize$gsize[match(pathway,
                table = gsSize$pathway
            )],
            propGenes = (lesize / gsize) * 100
        )

    return(plotDF)
}

# %%
read_gmt <- function(path) {
  con <- file(path, open = "r"); on.exit(close(con))
  parts <- strsplit(readLines(con, warn = FALSE), "\t", fixed = TRUE)
  sets <- lapply(parts, function(x) {
    genes <- x[-c(1, 2)]
    genes <- genes[genes != ""]
    unique(genes)
  })
  names(sets) <- toupper(vapply(parts, function(x) x[[1]], character(1)))
  sets
}

gmt_files <- c(
  "../../../data/gmt/MSigDB_Hallmark_2020.gmt",
  "../../../data/gmt/Reactome_Pathways_2024.gmt",
  "../../../data/gmt/ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.gmt",
  "../../../data/gmt/KEGG_2021_Human.gmt"
)
gmt_lists <- setNames(lapply(gmt_files, read_gmt), basename(gmt_files))

# %%
one_gmt_gsea_run <- function(gmx_list, label, cores = 8) {
  res <- RunGSEACelltype(
    logfc_list = deg_result,
    gmx        = gmx_list,
    ct.col     = "celltype",
    rank.col   = "rank_metric",
    gene.col   = "gene",
    collapsePathways = FALSE,
    ncores     = cores
  )
  dplyr::mutate(res, gmt_db = label)
}

# %% jupyter={"outputs_hidden": true}
dir.create("../../../data/rna/pseudobulk/results/fgsea_walds_results", recursive = TRUE, showWarnings = FALSE)
dir.create("../../../data/rna/pseudobulk/results/fgsea_walds_pvals", recursive = TRUE, showWarnings = FALSE)

files <- Sys.glob("../../../data/rna/pseudobulk/results/deseq2_results/deseq2*.csv")
# files <- files[!grepl("flu|bmmc", files)]

for (f in files) {
  message("Processing: ", f)
  base <- basename(f)
  stem <- sub("\\.csv$", "", base)
  stem <- sub("results_", "", stem)
  fgsea_csv <- sub("^results/deseq2_results/deseq2_results_", "../../../data/rna/pseudobulk/results/fgsea_walds_results/fgsea_walds_results_", f)
  pval_png <- sub("^results/deseq2_results/deseq2_results_", "../../../data/rna/pseudobulk/results/fgsea_walds_pvals/fgsea_walds_results_", sub("\\.csv$", "_pvals.png", f))

  # Extract condition for plot title
  condition <- sub(".*deseq2_results_(bmmc_|pbmc_)?(.*)\\.csv$", "\\2", f)

  # load DEGs
  deg_result <- data.table::fread(f)
  if (is.null(deg_result) || nrow(deg_result) == 0) {
    message("Skipping empty or null table: ", f)
    next
  }

  # p-value histogram per celltype
  options(repr.plot.width = 10, repr.plot.height = 6)
  p <- deg_result %>%
    ggplot(aes(x = pvalue, fill = celltype)) +
    geom_histogram(bins = 50) +
    facet_wrap(vars(celltype), scales = "free_y") +
    theme(legend.position = "none") +
    ggtitle(condition)
  ggsave(pval_png, plot = p, width = 16, height = 12)

  # rank metric -- change to walds
  deg_result <- deg_result %>%
      mutate(rank_metric = ifelse(is.na(padj), NA, stat)) %>%
      filter(!is.na(rank_metric))

  # run GSEA across all GMTs
  res_list <- purrr::imap(gmt_lists, ~ one_gmt_gsea_run(.x, .y, cores = 8))
  gsea_all <- dplyr::bind_rows(res_list)

  # save (swap _deg -> _fgsea)
  data.table::fwrite(gsea_all, fgsea_csv)

  message("Saved: ", pval_png, " and ", fgsea_csv)
}
