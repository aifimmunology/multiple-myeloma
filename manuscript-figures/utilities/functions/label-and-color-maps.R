# Author: Aishwarya Chander
# Organization: Allen Institute for Immunology

# ==============================================================================
# Centralized Color Maps and Label Definitions
# Purpose: Single source of truth for all colors and mappings used in figures
# Usage: source("01-scripts/label-and-color-maps.R") at the top of your script
# ==============================================================================

# ==============================================================================
# TIMEPOINT COLORS
# ==============================================================================

tp_color_map <- c(
  "PreTx"    = "#8B0000",
  "PI2C"     = "#C20272",
  "EI"       = "#7502C2",
  "ASCT60d"  = "#0021B0",
  "ASCT90d"  = "#2346DE",
  "ASCT1y"   = "#038AFF",
  "ASCT2y"   = "#008080",
  "Healthy"  = "#00470D"
)

# ==============================================================================
# SUBJECT COLORS
# ==============================================================================

subject_color_map <- c(
  "FH1001" = "#1f77b4",
  "FH1002" = "#ff7f0e",
  "FH1003" = "#279e68",
  "FH1004" = "#d62728",
  "FH1005" = "#aa40fc",
  "FH1006" = "#8c564b",
  "FH1007" = "#e377c2",
  "FH1008" = "#b5bd61",
  "FH1009" = "#17becf",
  "FH1010" = "#aec7e8",
  "FH1011" = "#ffbb78",
  "FH1012" = "#98df8a",
  # 'FH1013' = '',   # MISSING
  "FH1014" = "#ff9896",
  # 'FH1015' = '',   # MISSING
  "FH1016" = "#c5b0d5",
  "FH1017" = "#c49c94",
  "FH1018" = "#f7b6d2",
  "FH1019" = "#32CD32",
  "FH1020" = "#FF1493",
  "FH1021" = "#dbdb8d",
  "FH1022" = "#FFD700",
  "FH1023" = "#8B008B",
  "FH1024" = "#20B2AA",
  # 'FH1025' = '',   # MISSING
  "FH1026" = "#FF6347",
  "FH1027" = "#FF8C00",
  "FH1028" = "#BA55D3"
)

# ==============================================================================
# CELL TYPE LEVEL 3 COLORS
# ==============================================================================

l3_color_map <- c(
  # PBMC CELL TYPES
  # CD4 T cells
  "CD4 T Naive Core" = "#597FC6",
  "CD4 T Naive ISG+" = "#072E77",
  "CD4 T Naive SOX4+" = "#2B5777",
  "CD4 T CM" = "#6793A4",
  "CD4 T EM GZMB- CD27+" = "#406B9C",
  "CD4 T EM GZMB- CD27-" = "#233C5B",
  "CD4 T Mem ISG+" = "#5480A3",
  "CD4 T Mem KLRF1- GZMB+" = "#628088",

  # CD4 Tregs
  "Treg Naive" = "#8EB5D7",
  "Treg CD4 Mem" = "#001A9A",
  "Treg CD4 Mem KLRB1+" = "#3B76C1",
  "Treg CD4 Mem GZMK+" = "#A3B7C2",
  "Treg CD8 Mem" = "#C3D1DA",
  "Treg CD8 Mem KLRB1+" = "#C7E0C6",

  # CD8 T cells
  "CD8 T Naive Core" = "#63A686",
  "CD8 T Naive ISG+" = "#869A94",
  "CD8 T Naive SOX4+" = "#465F32",
  "CD8 T CM" = "#6D9F5E",
  "CD8 T EM GZMK+ CD27+" = "#95BDA1",
  "CD8 T EM KLRF1- GZMB+" = "#59895F",
  "CD8 T EM GZMK- CD27+" = "#314937",
  "CD8 T EM KLRF1+ GZMB+" = "#4E695B",
  "CD8 T Mem ISG+" = "#94C0AD",
  "CD8aa" = "#7ECDAD",

  # Other T cells
  "DN T" = "#656E68",
  "Prolif T" = "#AEAEB0",

  # gdT cells
  "gdT Vd1 SOX4+" = "#B6D7C8",
  "gdT Vd1 Naive" = "#537E41",
  "gdT Vd1 Eff KLRF1-" = "#A39E78",
  "gdT Vd1 Eff KLRF1+" = "#7CB38E",
  "gdT Vd2 GZMK+" = "#72AF5A",
  "gdT Vd2 GZMB+" = "#AECFB1",

  # MAIT cells
  "CD8 MAIT" = "#1E8F64",
  "CD4 MAIT" = "#6F723F",
  "MAIT ISG+" = "#BBB788",

  # B cells - PBMC
  "Trans B" = "#9C6469",
  "Naive B Core" = "#D5ABAB",
  "Naive B ISG+" = "#D9A996",
  "Mem B Early" = "#362027",
  "Mem B Core" = "#B16B72",
  "Mem B Activated" = "#60545A",
  "Mem B Type2" = "#4E2927",
  "Mem B CD95" = "#EBC3A7",
  "Eff B CD27+" = "#F4CEC3",
  "Eff B CD27-" = "#A3787B",
  "Plasma" = "#784D47",

  # NK cells
  "CD56br NK" = "#FF7F55",
  "CD56dim NK GZMK+" = "#A31001",
  "CD56dim NK GZMK-" = "#E16040",
  "NK Adaptive" = "#781B15",
  "CD56dim NK ISG+" = "#EEA172",
  "Prolif NK" = "#CC3B18",

  # Monocytes
  "CD14 Mono Core" = "#F6831B",
  "CD14 Mono ISG+" = "#DE5A07",
  "CD14 Mono IL1B+" = "#BE4C26",
  "Int Mono" = "#B75228",
  "CD16 Mono Core" = "#FCAF87",
  "CD16 Mono ISG+" = "#DB6A2E",
  "CD16 Mono C1Q+" = "#F4BA81",

  # Dendritic cells
  "cDC2 CD14+" = "#FEDC7A",
  "cDC2 HLA-DRhi" = "#BE6E23",
  "cDC2 ISG+" = "#8C6527",
  "cDC1" = "#F8C755",
  "pDC" = "#FBB64E",
  "ASDC" = "#8A4E1C",

  # Other PBMC
  "ILC" = "#F4AD3A",
  "CLP" = "#FCDB97",
  "CMP" = "#5B3006",
  "BaEoMaP" = "#3D1F15",
  "Platelet" = "#5B3930",

  # BMMC-SPECIFIC CELL TYPES
  # Pre B cells
  "Pre B Heavy" = "#bd71b5",
  "Pre B Light" = "#d08cc8",
  "Pre B ISG+" = "#a45c9d",
  "Pre B Prolif" = "#9a4892",

  # Prog B cells
  "Prog B" = "#C48AE7",
  "Prog B Pre" = "#b577d4",
  "Prog B Mature" = "#d39ff2",

  # Additional B cells from BMMC
  "Trans B Core" = "#9C6469",
  "Trans B ISG+" = "#B37980",

  # NK cells - additional BMMC
  "NK Tissue Res" = "#99362C",
  "NK Effector" = "#B3391F",

  # Monocytes - additional BMMC
  "Pre Mono Core" = "#B27FD9",
  "Pre Mono Prolif" = "#9D69C4",

  # Dendritic cells - additional BMMC
  "gdT" = "#8DC9A8",
  "cDC2 Core" = "#F0C95A",
  "Prog DC cDC" = "#8F4DC7",
  "Prog DC pDC" = "#A561D9",

  # Stem/Progenitor cells - BMMC
  "HSPC Stem" = "#6429A6",
  "HSPC Multi" = "#7538B8",
  "HSPC Prolif" = "#5420A0",
  "LMPP" = "#7A3DB8",
  "CMP Core" = "#5B3006",
  "CMP Gran" = "#A165D4",
  "MEP" = "#44177A",

  # Erythroid/Megakaryocyte
  "Pre Prog Ery" = "#361063",
  "Prog Ery" = "#1F1515",
  "Prog Ery Prolif" = "#1F083A",
  "Prog MK" = "#531E8E",

  # T cells
  "Treg" = "#4477C9",
  "MAIT" = "#3A8A5A",
  "CD8 T EM2" = "#2F5938",
  "CD8 T Tissue Res" = "#2D8069",
  "CD4 T EM1" = "#3F6BA0",
  "CD4 T EM2" = "#22507A",
  "CD4 T Mem Core" = "#5A8EA0",
  "CD8 T EM1" = "#537A50"
)

# ==============================================================================
# CELL TYPE LEVEL 3 TO LEVEL 1 MAPPING
# ==============================================================================

l3_to_l1 <- c(
  # --- Adaptive B ---
  "Eff B CD27+" = "Adaptive B",
  "Eff B CD27-" = "Adaptive B",
  "Mem B Activated" = "Adaptive B",
  "Mem B CD95" = "Adaptive B",
  "Mem B Core" = "Adaptive B",
  "Mem B Early" = "Adaptive B",
  "Mem B Type2" = "Adaptive B",

  # --- Adaptive Cytotoxic ---
  "CD56br NK" = "Adaptive Cytotoxic",
  "CD56dim NK GZMK+" = "Adaptive Cytotoxic",
  "CD56dim NK GZMK-" = "Adaptive Cytotoxic",
  "CD56dim NK ISG+" = "Adaptive Cytotoxic",
  "NK Adaptive" = "Adaptive Cytotoxic",
  "NK Effector" = "Adaptive Cytotoxic",
  "NK Tissue Res" = "Adaptive Cytotoxic",
  "CD8 T CM" = "Adaptive Cytotoxic",
  "CD8 T EM1" = "Adaptive Cytotoxic",
  "CD8 T EM2" = "Adaptive Cytotoxic",
  "CD8 T EM GZMK+ CD27+" = "Adaptive Cytotoxic",
  "CD8 T EM GZMK- CD27+" = "Adaptive Cytotoxic",
  "CD8 T EM KLRF1+ GZMB+" = "Adaptive Cytotoxic",
  "CD8 T EM KLRF1- GZMB+" = "Adaptive Cytotoxic",
  "CD8 T Mem ISG+" = "Adaptive Cytotoxic",
  "CD8 T Tissue Res" = "Adaptive Cytotoxic",
  "CD8aa" = "Adaptive Cytotoxic",
  "DN T" = "Adaptive Cytotoxic",
  "Prolif NK" = "Adaptive Cytotoxic",
  "Prolif T" = "Adaptive Cytotoxic",
  "gdT Vd1 Eff KLRF1+" = "Adaptive Cytotoxic",
  "gdT Vd1 Eff KLRF1-" = "Adaptive Cytotoxic",
  "gdT Vd1 Naive" = "Adaptive Cytotoxic",
  "gdT Vd1 SOX4+" = "Adaptive Cytotoxic",
  "gdT Vd2 GZMB+" = "Adaptive Cytotoxic",
  "gdT Vd2 GZMK+" = "Adaptive Cytotoxic",
  "gdT" = "Adaptive Cytotoxic",

  # --- Adaptive Helper ---
  "CD4 T CM" = "Adaptive Helper",
  "CD4 T EM1" = "Adaptive Helper",
  "CD4 T EM2" = "Adaptive Helper",
  "CD4 T EM GZMB- CD27+" = "Adaptive Helper",
  "CD4 T EM GZMB- CD27-" = "Adaptive Helper",
  "CD4 T Mem Core" = "Adaptive Helper",
  "CD4 T Mem ISG+" = "Adaptive Helper",
  "CD4 T Mem KLRF1- GZMB+" = "Adaptive Helper",

  # --- Immature Lymphoid ---
  "CD4 T Naive Core" = "Immature Lymphoid",
  "CD4 T Naive ISG+" = "Immature Lymphoid",
  "CD4 T Naive SOX4+" = "Immature Lymphoid",
  "CD8 T Naive Core" = "Immature Lymphoid",
  "CD8 T Naive ISG+" = "Immature Lymphoid",
  "CD8 T Naive SOX4+" = "Immature Lymphoid",
  "Naive B Core" = "Immature Lymphoid",
  "Naive B ISG+" = "Immature Lymphoid",
  "Trans B" = "Immature Lymphoid",
  "Trans B Core" = "Immature Lymphoid",
  "Trans B ISG+" = "Immature Lymphoid",
  "Treg" = "Immature Lymphoid",

  # --- Innate Lymphoid ---
  "MAIT" = "Innate Lymphoid",
  "CD4 MAIT" = "Innate Lymphoid",
  "CD8 MAIT" = "Innate Lymphoid",
  "ILC" = "Innate Lymphoid",
  "MAIT ISG+" = "Innate Lymphoid",

  # --- Mature Myeloid DC ---
  "cDC1" = "Mature Myeloid DC",
  "cDC2 Core" = "Mature Myeloid DC",
  "cDC2 CD14+" = "Mature Myeloid DC",
  "cDC2 HLA-DRhi" = "Mature Myeloid DC",
  "cDC2 ISG+" = "Mature Myeloid DC",
  "pDC" = "Mature Myeloid DC",
  "ASDC" = "Mature Myeloid DC",

  # --- Mature Myeloid Monocyte ---
  "CD14 Mono Core" = "Mature Myeloid Monocyte",
  "CD14 Mono IL1B+" = "Mature Myeloid Monocyte",
  "CD14 Mono ISG+" = "Mature Myeloid Monocyte",
  "CD16 Mono C1Q+" = "Mature Myeloid Monocyte",
  "CD16 Mono Core" = "Mature Myeloid Monocyte",
  "CD16 Mono ISG+" = "Mature Myeloid Monocyte",
  "Int Mono" = "Mature Myeloid Monocyte",

  # --- Mature Myeloid Platelet ---
  "Platelet" = "Mature Myeloid Platelet",

  # --- Plasma ---
  "Plasma" = "Plasma",

  # --- Progenitor B ---
  "Pre B Heavy" = "Progenitor B",
  "Pre B ISG+" = "Progenitor B",
  "Pre B Light" = "Progenitor B",
  "Pre B Prolif" = "Progenitor B",
  "Prog B" = "Progenitor B",
  "Prog B Mature" = "Progenitor B",
  "Prog B Pre" = "Progenitor B",

  # --- Progenitor DC ---
  "Prog DC cDC" = "Progenitor DC",
  "Prog DC pDC" = "Progenitor DC",

  # --- Progenitor HSPC ---
  "CLP" = "Progenitor HSPC",
  "LMPP" = "Progenitor HSPC",
  "HSPC Multi" = "Progenitor HSPC",
  "HSPC Prolif" = "Progenitor HSPC",
  "HSPC Stem" = "Progenitor HSPC",

  # --- Progenitor Monocyte ---
  "Pre Mono Core" = "Progenitor Monocyte",
  "Pre Mono Prolif" = "Progenitor Monocyte",

  # --- Progenitor Myelo-Ery ---
  "BaEoMaP" = "Progenitor Myelo-Ery",
  "CMP" = "Progenitor Myelo-Ery",
  "CMP Core" = "Progenitor Myelo-Ery",
  "CMP Gran" = "Progenitor Myelo-Ery",
  "MEP" = "Progenitor Myelo-Ery",
  "Pre Prog Ery" = "Progenitor Myelo-Ery",
  "Prog Ery" = "Progenitor Myelo-Ery",
  "Prog Ery Prolif" = "Progenitor Myelo-Ery",
  "Prog MK" = "Progenitor Myelo-Ery",

  # --- Regulatory T ---
  "Treg CD4 Mem" = "Regulatory T",
  "Treg CD4 Mem GZMK+" = "Regulatory T",
  "Treg CD4 Mem KLRB1+" = "Regulatory T",
  "Treg CD8 Mem" = "Regulatory T",
  "Treg CD8 Mem KLRB1+" = "Regulatory T",
  "Treg Naive" = "Regulatory T"
)

# ==============================================================================
# CELL TYPE LEVEL 1 COLORS
# ==============================================================================

l1_color_map <- c(
    # --- Lymphoid Lineage ---
    "Adaptive Helper" = "#0050A0",      # Deep Royal Blue (CD4 T)
    "Adaptive Cytotoxic" = "#008080",   # Deep Cyan/Teal (CD8 T, NK, gdT)
    "Immature Lymphoid" = "#75B0E0",    # Light Sky Blue (Naive/Transitional)
    "Innate Lymphoid" = "#50C878",      # Emerald Green (MAIT/ILC)
    "Regulatory T" = "#964B00",         # Darker Bronze/Brown 

    "Adaptive B" = "#D82243",           # Vivid Red (Memory/Effector B)
    "Plasma" = "#EA709A",               # Brighter Rose Pink
    
    # --- Myeloid Lineage ---
    "Mature Myeloid Monocyte" = "#FF7F00", # Vibrant Pumpkin Orange (Monocytes)
    "Mature Myeloid DC" = "#FFC300",    # Bright Gold
    "Mature Myeloid Platelet" = "#684C3F", # Rich Earth Brown (Platelets)

    # --- Progenitor Lineage (Improved Resolution) ---
    "Progenitor HSPC" = "#6A05A6",      # Deep Grape Purple (Stem Cells - Strong & distinct)
    "Progenitor B" = "#A070E0",         # Medium Lavender (B Progenitor - Distinct lighter purple)
    "Progenitor DC" = "#483D8B",        # Slate Blue-Purple (DC Progenitor - New blue-purple hue)
    "Progenitor Monocyte" = "#9A50BB",  # Orchid Purple (Monocyte Progenitor - Warmer, slightly pinker purple)
    "Progenitor Myelo-Ery" = "#E066FF"  # Vivid Magenta-Pink (Myelo-Erythroid Progenitor - Brightest, most distinct pink-purple)
)


l2_cmap <- c(
  # T cells
  "CD4 T Naive" = "#597FC6",
  "CD4 T Memory" = "#233C5B",
  "Treg" = "#001A9A",
  "DN T" = "#656E68",
  "Prolif T" = "#AEAEB0",
  "CD8 T Naive" = "#63A686",
  "CD8 T Memory" = "#314937",
  "CD8aa" = "#7ECDAD",
  "gdT" = "#72AF5A",
  "MAIT" = "#1E8F64",
  "Tissue Res T" = "#4A6B5C",

  # B cells
  "Trans B" = "#9C6469",
  "Naive B" = "#D5ABAB",
  "Memory B" = "#B16B72",
  "Effector B" = "#F4CEC3",
  "Plasma" = "#784D47",
  "Pre B" = "#bd71b5",
  "Prog B" = "#C48AE7",

  # NK cells
  "CD56br NK" = "#FF7F55",
  "CD56dim NK" = "#E16040",
  "Prolif NK" = "#CC3B18",
  "Tissue Res NK" = "#99362C",

  # Monocytes
  "CD14 Mono" = "#F6831B",
  "Int Mono" = "#B75228",
  "CD16 Mono" = "#FCAF87",
  "Pre Mono" = "#B27FD9",

  # Dendritic cells
  "cDC2" = "#BE6E23",
  "cDC1" = "#F8C755",
  "pDC" = "#FBB64E",
  "ASDC" = "#8A4E1C",
  "Prog DC" = "#8F4DC7",

  # Progenitors / Precursors
  "CLP" = "#A86CD9",
  "CMP" = "#9052C9",
  "LMPP" = "#7A3DB8",
  "HSPC" = "#6429A6",
  "Prog MK" = "#531E8E",
  "MEP" = "#44177A",
  "Pre Prog Ery" = "#361063",
  "Prog Ery" = "#2A0B4D",
  "Progenitor" = "#A31001",

  # Other
  "BaEoMaP" = "#A84C6A",
  "ILC" = "#F4AD3A",
  "Platelet" = "#5B3930"
)

# ==============================================================================
# FLOW CELL TYPE COLORS
# ==============================================================================
flow_color_map <- c(
  "CD14 Monocytes"                      = "#F6831B",
  "CD16 Monocytes"                      = "#FCAF87",
  "CD4 rm"                              = "#5480A3",
  "CD56 dim NK cell"                    = "#E16040",
  "CD56bright NK cell"                  = "#FF7F55",
  "CD8 rm"                              = "#94C0AD",
  "CM CD4+ T cell"                      = "#6793A4",
  "CM CD8 T cell"                       = "#6D9F5E",
  "Class-Switch Mem B"                  = "#60545A",
  "DN T cells"                          = "#656E68",
  "EM CD4 T cell"                       = "#406B9C",
  "EM CD8 T cell"                       = "#95BDA1",
  "Effector B"                          = "#F4CEC3",
  "Intermediate Monocytes"              = "#B75228",
  "Memory CD4 Treg"                     = "#001A9A",
  "Naive B"                             = "#D5ABAB",
  "Naive CD4 T cell"                    = "#597FC6",
  "Naive CD4 Treg"                      = "#8EB5D7",
  "Naive CD8 T cell"                    = "#63A686",
  "Non-Class Switched Mem B"            = "#B16B72",
  "Plasma cells"                        = "#784D47",
  "Plasmacytoid DC"                     = "#FBB64E",
  "Progenitor cell (circulating CD34+)" = "#6429A6",
  "T follicular helper"                 = "#2B5777",
  "TCR yd"                              = "#72AF5A",
  "TEMRA CD4"                           = "#233C5B",
  "TEMRA CD8 T cell"                    = "#59895F",
  "Transitional B cells"                = "#9C6469",
  "cDC1"                                = "#F8C755",
  "cDC2"                                = "#BE6E23"
)

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Get timepoint colors for a subset of visits
#' @param visits Character vector of visit names
#' @return Named vector of colors for the specified visits
get_timepoint_colors <- function(visits) {
  return(tp_color_map[visits])
}

#' Get subject colors for a subset of subjects
#' @param subjects Character vector of subject IDs
#' @return Named vector of colors for the specified subjects
get_subject_colors <- function(subjects) {
  return(subject_color_map[subjects])
}

#' Get L3 cell type colors for a subset of cell types
#' @param celltypes Character vector of cell type names
#' @return Named vector of colors for the specified cell types
get_celltype_colors <- function(celltypes) {
  return(l3_color_map[celltypes])
}

#' Get L1 category for L3 cell types
#' @param l3_celltypes Character vector of L3 cell type names
#' @return Character vector of L1 categories
get_l1_category <- function(l3_celltypes) {
  return(l3_to_l1[l3_celltypes])
}
