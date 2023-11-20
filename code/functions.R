# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# #### FUNCTIONS ###################

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CREATE DIRECTORY ####

#' Create a directory for project results
#'
#' @param data_dir The base data directory
#' @param project_name The name of the project
#' @param version The version code for the project
#' @param date The current date
#'
#' @return The path to the created directory
#'
create_directory <- function(result_dir = RESULT_DIR,
                                     project_name = SOFT_NAME,
                             version = VERSION_CODE,
                                     date = DATE) {
  dir_path <- paste0(result_dir, project_name, version, "_results_", date, "/")
  
  if (file.exists(dir_path)) {
    msg <- paste0("Directory '", dir_path, "' already existed. Check function's arguments")
  } else if (dir.create(dir_path)) {
    msg <- paste0("Directory '", dir_path, "' succesfully created.")
  } else {
    msg <- paste0("Directory ", dir_path ," cannot be created. Check function's arguments")
    stop(msg, call.=FALSE)
  } 
  message(msg)
  return(dir_path)
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CREATE DATA FIXED #############

#' Keeps observations that match the time condition 
#' 
#' @description
#' Selects observations for which there is information on the survival status at a fixed time instant, considering the relationship between survival status, follow-up time period, and the fixed time point.
#'
#' @param df_clinical A Clinical data frame object with the columns OS_STATUS and OS_MONTHS, and sample_id as row names.
#' @param int_time An integer defining the fixed time instant in days.
#'
#' @return A data frame object. The output has the following characteristics:
#' 
#' * Observations that match the time condition.
#' * Same columns as df_clinical with OS_STATUS and OS_MONTHS modified to the new time point.
#' 
#' @export
#'
#' @examples
create_data_fixed <- function(df_clinical, int_time) {
  tryCatch(
    expr = {
      # Check arguments
      if (!is.numeric(int_time)) {
        stop("Non-numeric time")
      } else if (int_time < 0) {
        stop("Non-positive time value")
      } else if (!is.data.frame(df_clinical)) {
        stop("Non-data frame object data")
      } else if (any(!(c("OS_MONTHS", "OS_STATUS") %in% colnames(df_clinical)))) {
        stop("Cannot find 'OS_MONTHS' or 'OS_STATUS' in column names")
      } else {
        ## Create a dataframe at the fixed time point
        clinical_data_fixed <- df_clinical
        
        ## Observations with a study time less than or equal to 1825 days and have not suffered the event; there is no information on the status at the fixed time, and the observation is eliminated.
        clinical_data_fixed <-
          clinical_data_fixed[-c(
            which(
              clinical_data_fixed$OS_MONTHS <= int_time &
                clinical_data_fixed$OS_STATUS == 0
            )
          ),]
        
        ## Observations with a study time of less than or equal to 1825 days and have suffered the event; retain this state at the fixed time.
        clinical_data_fixed$OS_STATUS[which(clinical_data_fixed$OS_MONTHS <= int_time &
                                              clinical_data_fixed$OS_STATUS == 1)] <- 1
        
        ## Observations with a study time greater than 1825 days and have suffered the event; they have not suffered it yet at the fixed date.
        clinical_data_fixed$OS_STATUS[which(clinical_data_fixed$OS_MONTHS > int_time &
                                              clinical_data_fixed$OS_STATUS == 1)] <- 0
        
        ## Observations with a study time greater than 1825 days and have not suffered the event; retain this state at the fixed int_time.
        clinical_data_fixed$OS_STATUS[which(clinical_data_fixed$OS_MONTHS > int_time &
                                              clinical_data_fixed$OS_STATUS == 0)] <- 0
        
        # Update the time of observations that exceed the set time
        clinical_data_fixed$OS_MONTHS[which(clinical_data_fixed$OS_MONTHS > int_time)] <-
          int_time
        
        ## Subset of the population followed up to a fixed time point
        clinical_data_fixed <-
          clinical_data_fixed[which(clinical_data_fixed$OS_MONTHS <= int_time),]
        
        return(clinical_data_fixed)
      }
    },
    error = function(e) {
      print(sprintf(
        "An error occurred in create_data_fixed at %s : %s",
        Sys.time(),
        e
      ))
    }
  )
  
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%
# MATCH SAMPLES ############

#' Select the intersection of observations in clinical and omic data frames.
#'
#' @param df_clinical Data frame: clinical data 
#' @param list_omics_df List of omics data frames
#'
#' @return A list of data frames with the same order as the input parameters. The set of observations is the same in each returned data frame.
#' @export
#'
#' @examples
match_samples <- function(df_clinical, list_omics_df) {
  tryCatch(
    expr = {
      # Check arguments
      if (!is.data.frame(df_clinical)) {
        error("Non-data frame type in 'df_clinical'")
      } else if (!is.list(list_omics_df)) {
        stop("Non-list type in 'list_omics_df'")
      } else {
        clinical_data <- df_clinical
        df_list <- list_omics_df
        dataframes_sample_ids <- list()
        dim_samples <- c()
        
        # Obtain sample_ids and dimensions of each dataframe
        for (i in 1:length(df_list)) {
          dataframes_sample_ids <- append(dataframes_sample_ids, list(colnames(df_list[i][[1]])))
          dim_samples <- append(dim_samples, length(colnames(df_list[i][[1]])))
        }
        
        # Since clinical data has sample_id as row names, we process it individually
        dataframes_sample_ids <- append(dataframes_sample_ids, list(row.names(clinical_data)))
        dim_samples <- append(dim_samples, length(row.names(clinical_data)))
        
        # Minimum vector of sample_ids
        samples_min <- dataframes_sample_ids[[which(dim_samples == min(dim_samples))[1]]]
        
        # Selecting from samples_min and all dataframes_samples_ids
        for (i in 1:length(dataframes_sample_ids)) {
          samples_min <- samples_min[which(samples_min %in% dataframes_sample_ids[i][[1]])]
        }
        
        # Update sample IDs
        for (i in 1:length(df_list)) {
          df_list[i][[1]] <- df_list[i][[1]][, c(which(colnames(df_list[i][[1]]) %in% samples_min))]
        }
        
        clinical_data <- clinical_data[c(which(row.names(clinical_data) %in% samples_min)), ]
        
        result <- list(clinical_data = clinical_data, omic_dfs = df_list)
        
        return(result)
      }
    },
    error = function(e) {
      print(sprintf("An error occurred in 'match_samples' at %s: %s",
                    Sys.time(),
                    e))
    }
  )
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# DIFFERENTIAL EXPRESSION ##########

#' Apply differential analysis to omics data frames
#'
#' @param df_clinical Data frame clinical_data derived of create_data_fixed
#' @param dfs_omics List of omics data frames
#' @param p_val Double type defining the P-value to be applied has filter and selecting the most representative variables
#'
#' @return A list of 3 lists. The characteristics of each list are the following:
#' 
#' * List 1 as ´diff_dfs´, store the omics data frames with the variables with p value < p_val
#' * List 2 as ´diff_names´, store a vector per each omic data frame of omic variables names which have passed the p value filter
#' * List 3 as ´top_table´, store a topTable element per each omic data frame with more information about it differential analysis 
#' associated
#' @export
#'
#' @examples
diff_express <- function(df_clinical, dfs_omics, p_val) {
  tryCatch(
    expr = {
      # Check arguments
      if (!is.data.frame(df_clinical)) {
        error("Non data.frame type in df_clinical")
      } else if (!is.list(dfs_omics)) {
        error("Non list type in df_clinical")
      } else if (!is.double(p_val)) {
        error("Non double type in p_val")
      } else {
        results <-
          list(
            diff_dfs = list(),
            diff_names = list(),
            top_table = list()
          )
        
        for (i in 1:length(dfs_omics)) {
          # Transform omics dataframes as data matrix
          omics_matrix <- as.matrix(dfs_omics[i][[1]])
          
          # Expression set
          express_set <- Biobase::ExpressionSet(log2(omics_matrix + 1))
          
          # Update patient data in expression set
          Biobase::pData(express_set) <-
            df_clinical[c(Biobase::sampleNames(express_set)), ]
          
          # Model matrix
          mm <-
            model.matrix(~ 0 + OS_STATUS, data = Biobase::pData(express_set))
          
          # lmFit
          fit <- limma::lmFit(express_set, mm)
          
          # makeContrast
          contrast <-
            limma::makeContrasts(OS_STATUS0 - OS_STATUS1, levels = colnames(coef(fit)))
          
          # contrast.fit
          contrast_fit <- limma::contrasts.fit(fit, contrast)
          
          # eBayes
          bayes <- limma::eBayes(contrast_fit)
          
          # topTable
          top_table <- limma::topTable(bayes, sort.by = "P", n = Inf, )
          
          # Filter topTable based on p value
          diff_names <-
            rownames(top_table[which(top_table$P.Value < p_val), ])
          
          # Update results
          results$diff_dfs[i][[1]] <-
            dfs_omics[i][[1]][c(diff_names), ]
          results$diff_names[i][[1]] <- diff_names
          results$top_table[i][[1]] <- top_table
        }
        
        return(results)
      }
    },
    error = function(e) {
      print(sprintf("An error occurred in diff_express at %s : %s",
                    Sys.time(),
                    e))
    }
  )
  
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# INFORMATION GAIN FILTER #######

#' Apply the information gain filter
#'
#' @param dataframe Data frame to which the filter is applied
#' @param int_n_features Integer defining the number of most representative features to be selected in the returned data frame
#'
#' @return A list of objects. Each object has the following characteristics:
#' 
#' * Object 1: df_ig, storing a data frame with the set of most representative features defined by 'int_n_features' as row names and sample_ids as column names.
#' * Object 2: result_ig, storing the information gain value associated with each variable in the original data frame 'dataframe'.
# @export
#'
#' @examples
informationGain_fs <- function(dataframe, int_n_features) {
  tryCatch(
    expr = {
      # Check arguments
      if (!is.data.frame(dataframe)) {
        stop("Argument 'dataframe' must be a data frame")
      } else if (!is.numeric(int_n_features)) {
        stop("Argument 'int_n_features' must be an integer")
      } else {
        ig_coreLearn <-
          attrEval(OS_STATUS == 1 ~ ., data = dataframe, estimator = "InfGain")
        
        df_result <-
          dataframe[, names(ig_coreLearn[order(ig_coreLearn, decreasing = TRUE)[1:int_n_features]])]
        df_result$OS_STATUS <- dataframe$OS_STATUS
        
        ig_return <-
          data.frame(Columns = names(ig_coreLearn[order(ig_coreLearn, decreasing = TRUE)]), IG = ig_coreLearn[order(ig_coreLearn, decreasing = TRUE)])
        
        elements <- list(df_ig = df_result, result_ig = ig_return)
        
        return(elements)
      }
    },
    error = function(e) {
      print(sprintf("An error occurred in informationGain_fs at %s : %s",
                    Sys.time(),
                    e))
    }
  )
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%
# EARLY INTEGRATION ##########


#' Apply early integration to a set of data frames
#'
#' @param list_dfs A list of data frames
#'
#' @return A data frame formed by concatenating variables from the datasets provided in 'list_dfs'.  
#' 
#' @export
#'
#' @examples
early_integration <- function(list_dfs) {
  tryCatch(
    expr = {
      # Check if list_dfs is a list of data frames
      if (!is.list(list_dfs) || !all(sapply(list_dfs, is.data.frame))) {
        stop("Argument 'list_dfs' must be a list of data frames")
      } else {
        dataframe <- list_dfs[[1]]
        
        for (i in 2:length(list_dfs)) {
          dataframe <- cbind(dataframe, list_dfs[[i]])
        }
        
        dataframe <- dataframe[!duplicated(as.list(dataframe))]
        
        return(dataframe)
      }
    },
    error = function(e) {
      print(sprintf("An error occurred in 'early_integration' at %s: %s",
                    Sys.time(),
                    e))
    }
  )
}

# %%%%%%%%%%%%%%%%%%%%%
# TRAIN MODEL #########

#' Training machine learning models applying hyperparameter grid search optimization
#'
#' @param char_model Character defining the model name
#' @param param_exp_grid expand.grid() element defining the hyperparameter space search 
#' @param df_data Data frame with the data 
#' @param list_cv_index List of train data indexes
#'
#' @return A list of objects. The characteristics of each object are as follows:
#' 
#' * Object 1: 'df_results', stores a data frame with model evaluation and bias prediction detection metrics results for each fold.
#' * Object 2: 'list_best_params', stores a list with the hyperparameters with the best results in each fold.
#' @export
#'
#' @examples
train_model <- function(char_model, param_exp_grid, df_data, list_cv_index) {
  tryCatch(
    expr = {
      if (!is.character(char_model)) {
        stop("Argument 'char_model' must be a character")
      } else if (!is.data.frame(df_data)) {
        stop("Argument 'df_data' must be a data frame")
      } else if (!is.list(list_cv_index) || !all(sapply(list_cv_index, is.numeric))) {
        stop("Argument 'list_cv_index' must be a list of numeric vectors")
      }
      
      df_results <- data.frame(Model = character(), Fold = numeric(), ACC = numeric(), AUC = numeric(), F1_score = numeric(), PRO_S = numeric(), PRO_NS = numeric(), PROPP_S = numeric(), PROPP_NS = numeric(), SEN_S = numeric(), SEN_NS = numeric(), EO_S = numeric(), EO_NS = numeric(), PRE_S = numeric(), PRE_NS = numeric(), PRP_S = numeric(), PRP_NS = numeric(), ACC_S = numeric(), ACC_NS = numeric(), ACCP_S = numeric(), ACCP_NS = numeric(), FNR_S = numeric(), FNR_NS = numeric(), FNRP_S = numeric(), FNRP_NS = numeric(), FPR_S = numeric(), FPR_NS = numeric(), FPRP_S = numeric(), FPRP_NS = numeric(), GZ_S = numeric(), GZ_NS = numeric())
      
      list_best_params <- list()
      
      for (i in 1:length(list_cv_index)) {
        training_data <- df_data[list_cv_index[[i]], ]
        levels(training_data$OS_STATUS) <- make.names(levels(factor(training_data$OS_STATUS)))
        test_data <- df_data[-list_cv_index[[i]], ]
        levels(test_data$OS_STATUS) <- make.names(levels(factor(test_data$OS_STATUS)))
        
        set.seed(3000)
        training_train_control <- trainControl(method = "cv", number = 2, search = "grid", classProbs = TRUE, summaryFunction = twoClassSummary)
        
        if (char_model == "svmLinear2") {
          training_model <- train(OS_STATUS ~ ., data = training_data, method = "svmLinear2", trControl = training_train_control, tuneGrid = param_exp_grid, metric = "ROC")
          
          final_model <- train(OS_STATUS ~ ., data = training_data, method = "svmLinear2", trControl = trainControl(method = "none", classProbs = TRUE), tuneGrid = training_model$bestTune, metric = "ROC")
        }
        
        if (char_model == "svmRadial") {
          training_model <- train(OS_STATUS ~ ., data = training_data, method = "svmRadialSigma", trControl = training_train_control, tuneGrid = param_exp_grid, metric = "ROC")
          
          final_model <- train(OS_STATUS ~ ., data = training_data, method = "svmRadialSigma", trControl = trainControl(method = "none", classProbs = TRUE), tuneGrid = training_model$bestTune, metric = "ROC")
        }
        
        if (char_model == "svmPoly") {
          training_model <- train(OS_STATUS ~ ., data = training_data, method = "svmPoly", trControl = training_train_control, tuneGrid = param_exp_grid, metric = "ROC")
          
          final_model <- train(OS_STATUS ~ ., data = training_data, method = "svmPoly", trControl = trainControl(method = "none", classProbs = TRUE), tuneGrid = training_model$bestTune, metric = "ROC")
        }
        
        test_predictions <- predict(final_model, test_data %>% select(-c("OS_STATUS")))
        
        cf_matrix_test <- confusionMatrix(test_predictions, test_data$OS_STATUS)
        accuracy_test <- as.vector(cf_matrix_test$overall[1])
        F1_test <- as.vector(cf_matrix_test$byClass[7])
        
        ## Binarize the outcome
        test_predictions <- as.character(test_predictions)
        test_predictions[test_predictions == "X0"] <- 0
        test_predictions[test_predictions == "X1"] <- 1
        test_data$OS_STATUS <- as.character(test_data$OS_STATUS)
        test_data$OS_STATUS[test_data$OS_STATUS == "X0"] <- 0
        test_data$OS_STATUS[test_data$OS_STATUS == "X1"] <- 1
        
        df_prediction <- data.frame(predictions = as.numeric(test_predictions), OS_STATUS = as.numeric(test_data$OS_STATUS), RACE = clinical_data_fixed$RACE[c(which(row.names(clinical_data_fixed) %in% row.names(test_data)))])
        # data(df_prediction)
        pred_AUC_test <- prediction(as.numeric(df_prediction$predictions), as.numeric(df_prediction$OS_STATUS))
        perf_AUC_test <- performance(pred_AUC_test, "auc")
        AUC_value_test <- perf_AUC_test@y.values[[1]]
        
        eq_odds <- equal_odds(
          data = df_prediction,
          outcome = "OS_STATUS",
          group = "RACE",
          preds = "predictions",
          base = "Black or African American"
        )
        
        eq_odds_sens_val <- eq_odds$Metric[1, 1]
        eq_odds_sens_par_val <- eq_odds$Metric[2, 1]
        eq_odds_nsens_val <- eq_odds$Metric[1, 2]
        eq_odds_nsens_par_val <- eq_odds$Metric[2, 2]
        
        pred_rt_parity <- pred_rate_parity(
          data = df_prediction,
          outcome = "OS_STATUS",
          group = "RACE",
          preds = "predictions",
          base = "Black or African American"
        )
        
        pred_rt_sens_val <- pred_rt_parity$Metric[1, 1]
        pred_rt_sens_par_val <- pred_rt_parity$Metric[2, 1]
        pred_rt_nsens_val <- pred_rt_parity$Metric[1, 2]
        pred_rt_nsens_par_val <- pred_rt_parity$Metric[2, 2]
        
        acc_par <- acc_parity(
          data = df_prediction,
          outcome = "OS_STATUS",
          group = "RACE",
          preds = "predictions",
          base = "Black or African American"
        )
        
        acc_sens_val <- acc_par$Metric[1, 1]
        acc_sens_par_val <- acc_par$Metric[2, 1]
        acc_nsens_val <- acc_par$Metric[1, 2]
        acc_nsens_par_val <- acc_par$Metric[2, 2]
        
        fnr_par <- fnr_parity(
          data = df_prediction,
          outcome = "OS_STATUS",
          group = "RACE",
          preds = "predictions",
          base = "Black or African American"
        )
        
        fnr_sens_val <- fnr_par$Metric[1, 1]
        fnr_sens_par_val <- fnr_par$Metric[2, 1]
        fnr_nsens_val <- fnr_par$Metric[1, 2]
        fnr_nsens_par_val <- fnr_par$Metric[2, 2]
        
        fpr_par <- fpr_parity(
          data = df_prediction,
          outcome = "OS_STATUS",
          group = "RACE",
          preds = "predictions",
          base = "Black or African American"
        )
        
        fpr_sens_val <- fpr_par$Metric[1, 1]
        fpr_sens_par_val <- fpr_par$Metric[2, 1]
        fpr_nsens_val <- fpr_par$Metric[1, 2]
        fpr_nsens_par_val <- fpr_par$Metric[2, 2]
        
        res_prop <- prop_parity(
          data = df_prediction,
          outcome = "OS_STATUS",
          group = "RACE",
          preds = "predictions",
          base = "Black or African American"
        )
        
        res_prop_sens_val <- res_prop$Metric[1, 1]
        res_prop_sens_par_val <- res_prop$Metric[2, 1]
        res_prop_nsens_val <- res_prop$Metric[1, 2]
        res_prop_nsens_par_val <- res_prop$Metric[2, 2]
        
        gz_s <- fpr_par$Metric[3, 1]
        gz_ns <- fpr_par$Metric[3, 2]
        
        
        df_results_add <- data.frame(
          Model = as.character(char_model),
          Fold = as.numeric(i),
          ACC = as.numeric(accuracy_test),
          AUC = as.numeric(AUC_value_test),
          F1_score = as.numeric(F1_test),
          PRO_S = as.numeric(res_prop_sens_val),
          PRO_NS = as.numeric(res_prop_nsens_val),
          PROPP_S = as.numeric(res_prop_sens_par_val),
          PROPP_NS = as.numeric(res_prop_nsens_par_val),
          SEN_S = as.numeric(eq_odds_sens_par_val),
          SEN_NS = as.numeric(eq_odds_nsens_par_val),
          EO_S = as.numeric(eq_odds_sens_par_val),
          EO_NS = as.numeric(eq_odds_nsens_par_val),
          PRE_S = as.numeric(pred_rt_sens_val),
          PRE_NS = as.numeric(pred_rt_nsens_val),
          PRP_S = as.numeric(pred_rt_sens_par_val),
          PRP_NS = as.numeric(pred_rt_nsens_par_val),
          ACC_S = as.numeric(acc_sens_val),
          ACC_NS = as.numeric(acc_nsens_val),
          ACCP_S = as.numeric(acc_sens_par_val),
          ACCP_NS = as.numeric(acc_nsens_par_val),
          FNR_S = as.numeric(fnr_sens_val),
          FNR_NS = as.numeric(fnr_nsens_val),
          FNRP_S = as.numeric(fnr_sens_par_val),
          FNRP_NS = as.numeric(fnr_nsens_par_val),
          FPR_S = as.numeric(fpr_sens_val),
          FPR_NS = as.numeric(fpr_nsens_val),
          FPRP_S = as.numeric(fpr_sens_par_val),
          FPRP_NS = as.numeric(fpr_nsens_par_val),
          GZ_S = as.numeric(gz_s),
          GZ_NS = as.numeric(gz_ns)
        )
        
        df_results <- rbind(df_results, df_results_add)
        
        list_best_params[i][[1]] <- training_model$bestTune
      }
      
      df_results_add_final <- data.frame(
        Model = as.character("Global"),
        Fold = NA,
        ACC = as.numeric(mean(df_results$ACC)),
        AUC = as.numeric(mean(df_results$AUC)),
        F1_score = as.numeric(mean(df_results$F1_score)),
        PRO_S = as.numeric(mean(df_results$PRO_S)),
        PRO_NS = as.numeric(mean(df_results$PRO_NS)),
        PROPP_S = as.numeric(mean(df_results$PROPP_S)),
        PROPP_NS = as.numeric(mean(df_results$PROPP_NS)),
        SEN_S = as.numeric(mean(df_results$SEN_S)),
        SEN_NS = as.numeric(mean(df_results$SEN_NS)),
        EO_S = as.numeric(mean(df_results$EO_S)),
        EO_NS = as.numeric(mean(df_results$EO_NS)),
        PRE_S = as.numeric(mean(df_results$PRE_S)),
        PRE_NS = as.numeric(mean(df_results$PRE_NS)),
        PRP_S = as.numeric(mean(df_results$PRP_S)),
        PRP_NS = as.numeric(mean(df_results$PRP_NS)),
        ACC_S = as.numeric(mean(df_results$ACC_S)),
        ACC_NS = as.numeric(mean(df_results$ACC_NS)),
        ACCP_S = as.numeric(mean(df_results$ACCP_S)),
        ACCP_NS = as.numeric(mean(df_results$ACCP_NS)),
        FNR_S = as.numeric(mean(df_results$FNR_S)),
        FNR_NS = as.numeric(mean(df_results$FNR_NS)),
        FNRP_S = as.numeric(mean(df_results$FNRP_S)),
        FNRP_NS = as.numeric(mean(df_results$FNRP_NS)),
        FPR_S = as.numeric(mean(df_results$FPR_S)),
        FPR_NS = as.numeric(mean(df_results$FPR_NS)),
        FPRP_S = as.numeric(mean(df_results$FPRP_S)),
        FPRP_NS = as.numeric(mean(df_results$FPRP_NS)),
        GZ_S = NA,
        GZ_NS = NA
      )
      
      df_results <- rbind(df_results, df_results_add_final)
      
      results <- list(df_results = df_results, list_best_params = list_best_params)
      
      return(results)
    },
    error = function(e) {
      print(sprintf("An error occurred in train_model at %s : %s",
                    Sys.time(),
                    e))
    }
  )
}

#' # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' # READ LGBM RESULTS INTO A LIST ########
#' 
#' #' Read LGBM results from files into a list
#' #'
#' #' @param base_path The base path for LGBM result files
#' #' @param num_files The number of LGBM result files to read
#' #' @param file_format The file format of LGBM result files (default is "csv")
#' #'
#' #' @return A list of data frames containing LGBM results
#' #' @export
#' #'
#' #' @examples
#' read_lgbm_results <- function(base_path, num_files, file_format = "csv") {
#'   tryCatch({
#'     # Check if base_path is a character string
#'     if (!is.character(base_path)) {
#'       stop("Argument 'base_path' must be a character string")
#'     }
#'     
#'     # Check if num_files is a numeric value
#'     if (!is.numeric(num_files) || length(num_files) != 1 || num_files <= 0) {
#'       stop("Argument 'num_files' must be a positive numeric value")
#'     }
#'     
#'     # Check if file_format is a character string
#'     if (!is.character(file_format) || !file_format %in% c("csv", "txt", "tsv")) {
#'       stop("Argument 'file_format' must be one of 'csv', 'txt', or 'tsv'")
#'     }
#'     
#'     # Function to read LGBM results given a base path and index
#'     read_lgbm_file <- function(index) {
#'       path <- paste0(base_path, "_", index, ".", file_format)
#'       read.delim(path, header = TRUE, sep = ifelse(file_format == "csv", ",", "\t"), dec = ".", fill = TRUE)
#'     }
#'     
#'     # Read LGBM results into a list
#'     list_lgbm_results <- lapply(0:(num_files - 1), read_lgbm_file)
#'     
#'     return(list_lgbm_results)
#'   }, error = function(e) {
#'     message(sprintf("An error occurred in read_lgbm_results at %s : %s", Sys.time(), e))
#'     return(NULL)
#'   })
#' }


# %%%%%%%%%%%%%%%%%%%%%
# BIAS LGBM ###########

#' Bias prediction detection applied to Light Gradient Boosting Machine (LGBM) Python results
#'
#' @param list_lgbm_results A data frame storing predictions and ground truth performed by LGBM
#' @param df_bias A data frame storing the sensible variable values (from clinical data fixed data frame)
#'
#' @return A data frame with model evaluation and bias prediction detection metrics results for each fold
#' @export
#'
#' @examples
bias_lgbm <- function(list_lgbm_results, df_bias) {
  
  tryCatch(
    expr = {
      if (!(all(sapply(list_lgbm_results, is.data.frame)))) {
        stop("Argument 'list_lgbm_results' must be a data frames list")
      } else {
        
        df_results <- data.frame(Model = as.character(), Fold = as.numeric(), PRO_S = as.numeric(), PRO_NS = as.numeric(), PROPP_S = as.numeric(), PROPP_NS = as.numeric(), SEN_S = as.numeric(), SEN_NS = as.numeric(), EO_S = as.numeric(), EO_NS = as.numeric(), PRE_S = as.numeric(), PRE_NS = as.numeric(), PRP_S = as.numeric(), PRP_NS = as.numeric(), ACC_S = as.numeric(), ACC_NS = as.numeric(), ACCP_S = as.numeric(), ACCP_NS = as.numeric(), FNR_S = as.numeric(), FNR_NS = as.numeric(), FNRP_S = as.numeric(), FNRP_NS = as.numeric(), FPR_S = as.numeric(), FPR_NS = as.numeric(), FPRP_S = as.numeric(), FPRP_NS = as.numeric(), GZ_S = as.numeric(), GZ_NS = as.numeric())
        
        
        for (i in 1:length(list_lgbm_results)) {
          list_lgbm_results[i][[1]] <- list_lgbm_results[i][[1]][c(which(list_lgbm_results[i][[1]]$X %in% row.names(df_bias))),]
          
          
          list_lgbm_results[i][[1]]$predictions <- factor(list_lgbm_results[i][[1]]$predictions)
          levels(list_lgbm_results[i][[1]]$predictions) <- list(X0 = 0, X1 = 1)
          list_lgbm_results[i][[1]]$truth <- factor(list_lgbm_results[i][[1]]$truth)
          levels(list_lgbm_results[i][[1]]$truth) <- list(X0 = 0, X1 = 1)
          
          df_prediction <- data.frame(predictions = list_lgbm_results[i][[1]]$predictions, OS_STATUS = list_lgbm_results[i][[1]]$truth, RACE = df_bias$RACE[c(which(row.names(df_bias) %in% list_lgbm_results[i][[1]]$X))])
          
          eq_odds <- equal_odds(
            data = df_prediction,
            outcome = "OS_STATUS",
            group = "RACE",
            preds = "predictions",
            base = "Black or African American"
          )
          
          eq_odds_sens_val <- eq_odds$Metric[1, 1]
          eq_odds_sens_par_val <- eq_odds$Metric[2, 1]
          eq_odds_nsens_val <- eq_odds$Metric[1, 2]
          eq_odds_nsens_par_val <- eq_odds$Metric[2, 2]
          
          pred_rt_parity <- pred_rate_parity(
            data = df_prediction,
            outcome = "OS_STATUS",
            group = "RACE",
            preds = "predictions",
            base = "Black or African American"
          )
          
          pred_rt_sens_val <- pred_rt_parity$Metric[1, 1]
          pred_rt_sens_par_val <- pred_rt_parity$Metric[2, 1]
          pred_rt_nsens_val <- pred_rt_parity$Metric[1, 2]
          pred_rt_nsens_par_val <- pred_rt_parity$Metric[2, 2]
          
          acc_par <- acc_parity(
            data = df_prediction,
            outcome = "OS_STATUS",
            group = "RACE",
            preds = "predictions",
            base = "Black or African American"
          )
          
          acc_sens_val <- acc_par$Metric[1, 1]
          acc_sens_par_val <- acc_par$Metric[2, 1]
          acc_nsens_val <- acc_par$Metric[1, 2]
          acc_nsens_par_val <- acc_par$Metric[2, 2]
          
          fnr_par <- fnr_parity(
            data = df_prediction,
            outcome = "OS_STATUS",
            group = "RACE",
            preds = "predictions",
            base = "Black or African American"
          )
          
          fnr_sens_val <- fnr_par$Metric[1, 1]
          fnr_sens_par_val <- fnr_par$Metric[2, 1]
          fnr_nsens_val <- fnr_par$Metric[1, 2]
          fnr_nsens_par_val <- fnr_par$Metric[2, 2]
          
          fpr_par <- fpr_parity(
            data = df_prediction,
            outcome = "OS_STATUS",
            group = "RACE",
            preds = "predictions",
            base = "Black or African American"
          )
          
          fpr_sens_val <- fpr_par$Metric[1, 1]
          fpr_sens_par_val <- fpr_par$Metric[2, 1]
          fpr_nsens_val <- fpr_par$Metric[1, 2]
          fpr_nsens_par_val <- fpr_par$Metric[2, 2]
          
          res_prop <- prop_parity(
            data = df_prediction,
            outcome = "OS_STATUS",
            group = "RACE",
            preds = "predictions",
            base = "Black or African American"
          )
          
          res_prop_sens_val <- res_prop$Metric[1, 1]
          res_prop_sens_par_val <- res_prop$Metric[2, 1]
          res_prop_nsens_val <- res_prop$Metric[1, 2]
          res_prop_nsens_par_val <- res_prop$Metric[2, 2]
          
          gz_s <- fpr_par$Metric[3, 1]
          gz_ns <- fpr_par$Metric[3, 2]
          
          
          df_results_add <- data.frame(
            Model = "LGBM",
            Fold = as.numeric(i),
            PRO_S = as.numeric(res_prop_sens_val),
            PRO_NS = as.numeric(res_prop_nsens_val),
            PROPP_S = as.numeric(res_prop_sens_par_val),
            PROPP_NS = as.numeric(res_prop_nsens_par_val),
            SEN_S = as.numeric(eq_odds_sens_par_val),
            SEN_NS = as.numeric(eq_odds_nsens_par_val),
            EO_S = as.numeric(eq_odds_sens_par_val),
            EO_NS = as.numeric(eq_odds_nsens_par_val),
            PRE_S = as.numeric(pred_rt_sens_val),
            PRE_NS = as.numeric(pred_rt_nsens_val),
            PRP_S = as.numeric(pred_rt_sens_par_val),
            PRP_NS = as.numeric(pred_rt_nsens_par_val),
            ACC_S = as.numeric(acc_sens_val),
            ACC_NS = as.numeric(acc_nsens_val),
            ACCP_S = as.numeric(acc_sens_par_val),
            ACCP_NS = as.numeric(acc_nsens_par_val),
            FNR_S = as.numeric(fnr_sens_val),
            FNR_NS = as.numeric(fnr_nsens_val),
            FNRP_S = as.numeric(fnr_sens_par_val),
            FNRP_NS = as.numeric(fnr_nsens_par_val),
            FPR_S = as.numeric(fpr_sens_val),
            FPR_NS = as.numeric(fpr_nsens_val),
            FPRP_S = as.numeric(fpr_sens_par_val),
            FPRP_NS = as.numeric(fpr_nsens_par_val),
            GZ_S = as.numeric(gz_s),
            GZ_NS = as.numeric(gz_ns)
          )
          
          df_results <- rbind(df_results, df_results_add)
        }
        
        # Result computations and return statement
        
        df_results_add_final <- data.frame(
          Model = as.character("Global"),
          Fold = NA,
          PRO_S = as.numeric(mean(df_results$PRO_S)),
          PRO_NS = as.numeric(mean(df_results$PRO_NS)),
          PROPP_S = as.numeric(mean(df_results$PROPP_S)),
          PROPP_NS = as.numeric(mean(df_results$PROPP_NS)),
          SEN_S = as.numeric(mean(df_results$SEN_S)),
          SEN_NS = as.numeric(mean(df_results$SEN_NS)),
          EO_S = as.numeric(mean(df_results$EO_S)),
          EO_NS = as.numeric(mean(df_results$EO_NS)),
          PRE_S = as.numeric(mean(df_results$PRE_S)),
          PRE_NS = as.numeric(mean(df_results$PRE_NS)),
          PRP_S = as.numeric(mean(df_results$PRP_S)),
          PRP_NS = as.numeric(mean(df_results$PRP_NS)),
          ACC_S = as.numeric(mean(df_results$ACC_S)),
          ACC_NS = as.numeric(mean(df_results$ACC_NS)),
          ACCP_S = as.numeric(mean(df_results$ACCP_S)),
          ACCP_NS = as.numeric(mean(df_results$ACCP_NS)),
          FNR_S = as.numeric(mean(df_results$FNR_S)),
          FNR_NS = as.numeric(mean(df_results$FNR_NS)),
          FNRP_S = as.numeric(mean(df_results$FNRP_S)),
          FNRP_NS = as.numeric(mean(df_results$FNRP_NS)),
          FPR_S = as.numeric(mean(df_results$FPR_S)),
          FPR_NS = as.numeric(mean(df_results$FPR_NS)),
          FPRP_S = as.numeric(mean(df_results$FPRP_S)),
          FPRP_NS = as.numeric(mean(df_results$FPRP_NS)),
          GZ_S = NA,
          GZ_NS = NA
        )
        
        df_results <- rbind(df_results, df_results_add_final)
        
        return(df_results)
      }
    },
    error = function(e) {
      print(sprintf("An error occurred in bias_lgbm at %s : %s", Sys.time(), e))
    }
  )
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# READ LGBM RESULTS INTO A LIST ########

#' Read LGBM results from files into a list
#'
#' @param base_path The base path for LGBM result files
#' @param num_files The number of LGBM result files to read
#' @param file_format The file format ("csv", "txt", "tsv") of LGBM result files (default is "csv")
#'
#' @return A list of data frames containing LGBM results
#' @export
#'
#' @examples
read_lgbm_results <- function(base_path, num_files, file_format = "csv") {
  tryCatch({
    # Check if base_path is a character string
    if (!is.character(base_path)) {
      stop("Argument 'base_path' must be a character string")
    }
    
    # Check if num_files is a numeric value
    if (!is.numeric(num_files) || length(num_files) != 1 || num_files <= 0) {
      stop("Argument 'num_files' must be a positive numeric value")
    }
    
    # Check if file_format is a character string
    if (!is.character(file_format) || !file_format %in% c("csv", "txt", "tsv")) {
      stop("Argument 'file_format' must be one of 'csv', 'txt', or 'tsv'")
    }
    
    # Construct the file path using the base path + "_" + index, and the file format
    read_lgbm_file <- function(index) {
      path <- paste0(base_path, "_", index, ".", file_format)
      read.delim(path, header = TRUE, sep = ifelse(file_format == "csv", ",", "\t"), dec = ".", fill = TRUE)
    }
    
    # Read LGBM results into a list
    list_lgbm_results <- lapply(0:(num_files - 1), read_lgbm_file)
    
    return(list_lgbm_results)
  }, error = function(e) {
    message(sprintf("An error occurred in read_lgbm_results at %s : %s", Sys.time(), e))
    return(NULL)
  })
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# TRANSLATE DATAEXPLORER COLUMNS NAMES ########

#' Translate column names of a DataExplorer::introduce data frame to Spanish
#'
#' This function takes a DataExplorer::introduce data frame with specific column names and translates them to
#' their equivalent names in Spanish, following a predefined mapping.
#'
#' Spanish output column names = "Filas", "Columnas", "Columnas Discretas", "Columnas Continuas", "Total Columnas Nulas", "NAs", "Observaciones Totales", "Filas Completas" 
#'
#' @param df_introduce The input DataExplorer::introduce data frame with original column names
#'
#' @return A data frame with column names translated to Spanish
#'
translate_introduce_columns <- function(df_introduce) {
  
  english_columns <- c("rows", "columns", "discrete_columns", "continuous_columns", "all_missing_columns", "total_missing_values", "complete_rows", "total_observations", "memory_usage")
  
  tryCatch({
    if (!(all(english_columns %in% colnames(df_introduce)))) {
      stop("Check input data frame. Not all the default DataExplorer::introduce data frame columns found.")
    } else {
      colnames(df_introduce)[which(names(df_introduce) == "rows")] <- "Filas"
      colnames(df_introduce)[which(names(df_introduce) == "columns")] <- "Columnas"
      colnames(df_introduce)[which(names(df_introduce) == "discrete_columns")] <- "Columnas Discretas"
      colnames(df_introduce)[which(names(df_introduce) == "continuous_columns")] <- "Columnas Continuas"
      colnames(df_introduce)[which(names(df_introduce) == "all_missing_columns")] <- "Total Columnas Nulas"
      colnames(df_introduce)[which(names(df_introduce) == "total_missing_values")] <- "NAs"
      colnames(df_introduce)[which(names(df_introduce) == "total_observations")] <- "Observaciones Totales"
      colnames(df_introduce)[which(names(df_introduce) == "complete_rows")] <- "Filas Completas"
      rm(english_columns)
    }
    
    return(df_introduce)
  }, error = function(e) {
    message(sprintf("An error occurred in translate_introduce_columns at %s : %s", Sys.time(), e))
    return(NULL)
  })
}

cat("\n", "############## FINISH LOADING FUNCTIONS ####################", "\n")