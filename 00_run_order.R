# Run order for the clean R code package
# Put this folder in the project root or update project_root paths inside scripts before running.
# Some analyses require processed data tables and sequence-processing outputs deposited separately.

# 1. QC and attrition summaries
# source('01_qc_and_attrition.R')

# 2. Effect-size tables and overall meta-analysis
# source('02_effect_size_and_overall_meta.R')

# 3. Categorical subgroup analyses and main forest plots
# source('03_categorical_subgroups_and_main_forest.R')

# 4. Continuous moderators and random-forest screening
# source('04_continuous_moderators_and_random_forest.R')

# 5. Dominant phylum composition analyses
# source('05_taxonomic_composition_phylum.R')

# 6. Network metric summaries and network figures
# source('06_network_analysis_and_figures.R')

# 7. Publication-bias diagnostics and sensitivity checks
# source('07_publication_bias_and_sensitivity.R')

# 8. Publication figure wrapper
# source('08_publication_figures.R')
