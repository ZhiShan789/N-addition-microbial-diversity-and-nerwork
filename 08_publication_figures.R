# Publication figure generation wrapper
# Clean GitHub code package for the N deposition soil microbiome synthesis.
# This file combines the R scripts actually used for the manuscript-level analyses.

# Sources combined: 1
# Generated on: 2026-06-25


# ================================================================
# Source: Data_analysis/Remote_working_package/06_R_scripts/run_all_publication_figures.R
# ================================================================

setwd("E:/BaiduSyncdisk/N_deposition1/Data_analysis")

source("analyze_phylum_composition_significance.R")
source("plot_network_paired_metrics_pretty.R")

message("All R-based publication figures have been regenerated.")

