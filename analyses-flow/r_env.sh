#!/usr/bin/bash

R_ENV="/home/workspace/environment/pbmc_flow_r_env"
# HISE-specific path — adjust manually if running locally

# Create conda environment with R and essential packages
conda create -y -p $R_ENV -c conda-forge \
    r-base=4.3 \
    r-ggplot2 r-dplyr r-tidyr r-ggalluvial \
    r-readr r-stringr r-lubridate \
    r-data.table r-magrittr \
    r-scales r-viridis r-rcolorbrewer r-patchwork r-cowplot r-forcats r-compositions \
    r-corrplot r-reshape2 \
    r-irkernel r-tidyverse \
    r-gridextra r-rlang \
    r-broom r-ggpubr

# Activate the environment
conda activate $R_ENV

# Install the R Jupyter kernel
R -e "IRkernel::installspec(name = 'pbmc_flow_r_env', displayname = 'R (PBMC Flow)')"