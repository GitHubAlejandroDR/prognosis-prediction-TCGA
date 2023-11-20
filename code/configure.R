# %%%%%%%%%%%%%%%%%%%%%%%
# ENVIRONMENT PREPARATION

rm(list=ls())
gc()
graphics.off()

PROJECT_NAME = "Applying machine learning to predict pathology response in colorectal cancer"

# %%%%%%%%%%%%%%%%%%%
# LIBRARIES ARGUMENTS

PKG_UPDATE = TRUE
VERBOSE_MODE = TRUE

# %%%%%%%%%%%%%%%%
# BASE DIRECTORIES

CODE_DIR = r"(C:/Users/aleja/Documents/github_repositories/prognosis-prediction-TCGA/code/)"
RMDs_DIR = r"(C:/Users/aleja/Documents/github_repositories/prognosis-prediction-TCGA/RMDs/)"
DATA_DIR = r"(C:/Users/aleja/Documents/github_repositories/prognosis-prediction-TCGA/data/raw/)"
DATA_PROCESSED_DIR <- r"(C:/Users/aleja/Documents/github_repositories/prognosis-prediction-TCGA/data/processed/)"
RESULT_DIR <- r"(C:/Users/aleja/Documents/github_repositories/prognosis-prediction-TCGA/results/)"
PLOTS_DIR <- r"(C:\Users\aleja\Documents\github_repositories\prognosis-prediction-TCGA\plots\)"

# %%%%%%%%%%%%%%
# RAW DATA PATHS

RAW_CLINICAL_PATIENT <- r"(C:/Users/aleja/Documents/github_repositories/prognosis-prediction-TCGA/data/raw/data_clinical_patient.txt)"
RAW_CLINICAL_SAMPLE <- r"(C:/Users/aleja/Documents/github_repositories/prognosis-prediction-TCGA/data/raw/data_clinical_sample.txt)"
RAW_GENOMIC <- r"(C:/Users/aleja/Documents/github_repositories/prognosis-prediction-TCGA/data/raw/data_mrna_seq_v2_rsem.txt)"
RAW_METHYLATION <- r"(C:/Users/aleja/Documents/github_repositories/prognosis-prediction-TCGA/data/raw/data_methylation.txt)"
RAW_MICROBIOME <- r"(C:/Users/aleja/Documents/github_repositories/prognosis-prediction-TCGA/data/raw/data_microbiome.txt)"

# %%%%%%%%%%%%%%%%%%%%%
# STATISTICAL CONSTANTS

P <- 0.05

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# MACHINE LEARNING CONFIGURATION PARAMETERS

PERC <- 20
FOLD <- 5
SEED <- 3000
NFEATURES <- 30
DAYS <- 1825


# ---------------------

cat("\n", "############## FINISH LOADING CONFIGURATION ####################", "\n")
