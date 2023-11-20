# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CHECKING EXECUTION SESSION #######

error_msg <- "### ATENTION - Current execution non interactive ####"

ARGS <- commandArgs(trailingOnly = TRUE)
if (length(ARGS) == 1) { 
  message("\n", "### Interactive execution. ARGS[1] will be trated as configuration file ###", "\n")
  source(ARGS[1])
} else if (length(ARGS) > 1) {
  warning(" ### ARGS in positions > 1 will be depreced ###")
} else if (!(interactive())) {
  stop(error_msg, call.=FALSE)
} 

# %%%%%%%%%%%%%%%%%%%%%%
# LOADING LIBRARIES ####


file_to_source <- paste0(CODE_DIR, "libraries.R")
source(file_to_source)


# %%%%%%%%%%%%%%%%%%%%%%
# LOADING FUNCTIONS ####


file_to_source <- paste0(CODE_DIR, "functions.R")
source(file_to_source)

rm(file_to_source)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CREATING EXECUTION RESULTS DIRECTORY #######

SHORT_NAME <- "TFG_"
VERSION_CODE <- 0.1
VERSION_R <- R.Version()
DATE <- format(Sys.time(), "%F_%H.%M.%S")

WD <- create_directory(RESULT_DIR, SHORT_NAME, VERSION_CODE, DATE)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CHECKING DATA_DIR ##########

if (!file.exists(DATA_DIR)) {
  msg <- paste0(" ### The directory ",
                DATA_DIR, " does not exist. Review the DATA_DIR variable ###\n")
  stop(msg, call. = FALSE)
}

setwd(WD)

cat("\n", "############## CREATING MARKDOWN REPORT ####################", "\n")

main_RMD <- paste0(CODE_DIR, "Report.Rmd")
render(main_RMD, output_dir = WD, quiet = TRUE)

