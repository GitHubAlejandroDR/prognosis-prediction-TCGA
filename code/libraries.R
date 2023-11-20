# %%%%%%%%%%%%%%%%%%%%%%%%%%
# ##### LIBRARIES ##########

# %%%%%%%%%%%%%%%%%%%%%
# CRAN libraries ######
libraries <- c("fairness", "dplyr", "CORElearn", "FSelector", "randomForest", "stringr", "caret", "splitTools", "MLmetrics", "ROCR", "glmnet", "kernlab", "mltools", "data.table", "AER", "tidyverse", "naniar", "gridExtra", "knitr", "rmarkdown", "markdown", "DataExplorer")

# %%%%%%%%%%%%%%%%%%%%%%%%%%
# Install CRAN libraries ###
installed_libraries <- libraries %in% rownames(installed.packages())
if (any(installed_libraries == FALSE)) {
  install.packages(libraries[!installed_libraries], force=T)
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%
# Load CRAN libraries ###########
invisible(lapply(libraries, library, character.only = TRUE))

# %%%%%%%%%%%%%%%%%%%%%%%%%%
# Install BiocManager #####
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# %%%%%%%%%%%%%%%%%%%%%%%%%%
# BiocManager libraries ####
bioc_libraries <- c("limma", "Biobase")

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Install BiocManager libraries ####
BiocManager::install(bioc_libraries, force = T)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load BiocManager libraries ####
invisible(lapply(bioc_libraries, library, character.only = TRUE))

cat("\n", "############## FINISH LOADING LIBRARIES ####################", "\n")

