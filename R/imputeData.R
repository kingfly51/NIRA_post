#' Advanced Dichotomous Data Imputation with Multiple Methods
#'
#' @description
#' Provides a unified interface for various statistical and machine learning imputation methods.
#' Supports both univariate and multivariate imputation with automatic method validation
#' and comprehensive error handling. Part of the DataPreprocessingR package ecosystem.
#'
#' @param data Data frame or matrix containing missing values to be imputed.
#' @param method Imputation method to use. One of:  
#'   `"mode"`, `"knn"`, `"random_forest"`, `"decision_tree"`, `"svm"`,  
#'   `"ann"`, `"mi"`, `"logistic"` (default: `"mi"`)
#' @param ... Additional arguments passed to underlying imputation functions
#'
#' @return
#' A complete dataset with missing values imputed, matching the input structure.
#' Returns the original data unchanged if no missing values are detected.
#'
#' @examples
#' \dontrun{
#' # Create dataset with missing values
#' df <- data.frame(
#'   A = c(1, 2, 3, 4, 5),
#'   B = c(NA, 2, NA, 4, 5),
#'   C = c(1, 2, 3, NA, NA)
#' )
#' 
#' # Check missing patterns first
#' missing_info <- check_missing(df)
#' 
#' # Perform imputation using random forest
#' imputed_data <- impute_data(df, method = "random_forest")
#' 
#' # Perform imputation with custom parameters
#' imputed_data <- impute_data(df, method = "knn", k = 3)
#' imputed_data <- impute_data(df, method = "mi", m = 10, maxit = 20)
#' }
#' 
#' @section Method Details:
#' - **`mode`**: Simple majority voting for categorical data
#' - **`knn`**: k-Nearest Neighbors using VIM package
#' - **`random_forest`**: MissForest algorithm for mixed data types
#' - **`decision_tree`**: CART-based imputation using rpart
#' - **`svm`**: Support Vector Machines for regression/classification
#' - **`ann`**: Artificial Neural Networks with nnet
#' - **`mi`**: Multiple Imputation by Chained Equations (MICE)
#' - **`logistic`**: Logistic regression for binary outcomes
#'
#' @note
#' The function automatically validates method suitability based on data completeness.
#' Methods requiring complete features (`knn`, `decision_tree`, `svm`, `ann`, `logistic`)
#' will fail if no complete columns exist. Use `check_missing()` to verify method availability.
#'
#' @author
#' Fei Wang, State Key Laboratory of Cognitive Science and Mental Health, 
#' Institute of Psychology, Chinese Academy of Sciences (2025-08-31)
#'
#' @section Requirements:
#' Feature-based methods (knn, decision_tree, svm, ann, logistic) require 
#' at least one complete column. Use `check_missing()` to diagnose data patterns
#' before imputation.
#'
#' @seealso
#' * [`check_missing()`] for data diagnostics before imputation
#' * [`mice::mice()`] for multiple imputation implementation
#' * [`missForest::missForest()`] for random forest imputation
#' * [`VIM::kNN()`] for k-nearest neighbors imputation
#' * [`rpart::rpart()`] for decision tree models
#' * [`e1071::svm()`] for support vector machines
#' * [`nnet::nnet()`] for neural networks implementation
#'
#' @importFrom stats complete.cases as.formula predict glm binomial
#' @importFrom VIM kNN
#' @importFrom missForest missForest
#' @importFrom mice mice complete
#' @importFrom rpart rpart
#' @importFrom e1071 svm
#' @importFrom nnet nnet
#' @export
imputeData <- function(data, method = "mi", ...) {
  # Input validation
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("Input must be a data frame or matrix")
  }
  
  if (sum(is.na(data)) == 0) {
    message("No missing values found in the data. Returning original data.")
    return(data)
  }
  
  # Check method validity
  valid_methods <- c("mode", "knn", "random_forest", "decision_tree", 
                     "svm", "ann", "mi", "logistic")
  
  if (!method %in% valid_methods) {
    stop("Invalid method. Available methods: ", paste(valid_methods, collapse = ", "))
  }
  
  # Check if method is appropriate for data pattern
  complete_cols <- names(which(colSums(is.na(data)) == 0))
  if (method %in% c("knn", "decision_tree", "svm", "ann", "logistic") && 
      length(complete_cols) == 0) {
    stop("Method '", method, "' requires at least one complete column. ",
         "Use check_missing() to see available methods.")
  }
  
  # Perform imputation based on selected method
  result <- switch(method,
                   "mode" = {
                     message("Using mode imputation...")
                     mode_data <- data
                     for (col in names(data)) {
                       if (any(is.na(data[[col]]))) {
                         mode_val <- names(sort(table(data[[col]]), decreasing = TRUE))[1]
                         imputed_data[is.na(data[[col]]), col] <- mode_val
                       }
                     }
                     mode_data
                   },
                   "knn" = {
                     message("Using k-Nearest Neighbors imputation (VIM::kNN)...")
                     dots <- list(...)
                     k <- ifelse("k" %in% names(dots), dots$k, 5)
                     knn_result <- VIM::kNN(data, k = k, ...)
                     imp_cols <- grep("_imp$", names(knn_result), value = TRUE)
                     if (length(imp_cols) > 0) {
                       knn_result <- knn_result[, !names(knn_result) %in% imp_cols]
                     }
                     knn_result <- knn_result[, names(data)]
                     knn_result
                   },
                   "random_forest" = {
                     message("Using Random Forest imputation...")
                     dots <- list(...)
                     maxiter <- ifelse("maxiter" %in% names(dots), dots$maxiter, 10)
                     ntree <- ifelse("ntree" %in% names(dots), dots$ntree, 100)
                     rf_result <- missForest::missForest(data, maxiter = maxiter, ntree = ntree, ...)
                     rf_result$ximp
                   },
                   "decision_tree" = {
                     message("Using Decision Tree imputation...")
                     imputed_data <- data
                     
                     for (col in names(which(colSums(is.na(data)) > 0))) {
                       message("  Imputing column: ", col)
                       
                       # Prepare training data
                       train_data <- data[complete.cases(data[[col]]), ]
                       train_data <- train_data[complete.cases(train_data[complete_cols]), ]
                       
                       if (nrow(train_data) <= 1) {
                         warning("No training data available for column: ", col, ". Skipping imputation.")
                         next
                       }
                       else if (nrow(train_data) > 1) {
                         # Build formula
                         formula <- as.formula(paste(col, "~", paste(complete_cols, collapse = "+")))
                         
                         # Train decision tree model
                         tree_model <- rpart::rpart(formula, data = train_data, method = "class")
                         
                         # Predict missing values
                         missing_indices <- which(is.na(data[[col]]))
                         predict_data <- data[missing_indices, complete_cols, drop = FALSE]
                         
                         # Ensure prediction data has no missing values
                         complete_predict <- complete.cases(predict_data)
                         if (any(complete_predict)) {
                           predictions <- predict(tree_model, 
                                                  predict_data[complete_predict, , drop = FALSE], 
                                                  type = "class")
                           imputed_data[missing_indices[complete_predict], col] <- as.numeric(as.character(predictions))
                         }
                       }
                     }
                     imputed_data
                   },
                   "svm" = {
                     message("Using SVM imputation...")
                     imputed_data <- data
                     
                     for (col in names(which(colSums(is.na(data)) > 0))) {
                       message("  Imputing column: ", col)
                       
                       # Prepare training data
                       train_data <- data[complete.cases(data[[col]]), ]
                       train_data <- train_data[complete.cases(train_data[complete_cols]), ]
                       
                       if (nrow(train_data) <= 1) {
                         warning("No training data available for column: ", col, ". Skipping imputation.")
                         next
                       }
                       else if (nrow(train_data) > 1) {
                         # Build formula
                         formula <- as.formula(paste(col, "~", paste(complete_cols, collapse = "+")))
                         
                         # Train SVM model
                         svm_model <- e1071::svm(formula, data = train_data)
                         
                         # Predict missing values
                         missing_indices <- which(is.na(data[[col]]))
                         predict_data <- data[missing_indices, complete_cols, drop = FALSE]
                         
                         # Ensure prediction data has no missing values
                         complete_predict <- complete.cases(predict_data)
                         if (any(complete_predict)) {
                           predictions <- predict(svm_model, 
                                                  predict_data[complete_predict, , drop = FALSE])
                           imputed_data[missing_indices[complete_predict], col] <- as.numeric(as.character(predictions))
                         }
                       }
                     }
                     imputed_data
                   },
                   "ann" = {
                     message("Using Neural Network imputation...")
                     dots <- list(...)
                     size <- ifelse("size" %in% names(dots), dots$size, 5)
                     maxit <- ifelse("maxit" %in% names(dots), dots$maxit, 200)
                     
                     imputed_data <- data
                     
                     for (col in names(which(colSums(is.na(data)) > 0))) {
                       message("  Imputing column: ", col)
                       
                       # Prepare training data
                       train_data <- data[complete.cases(data[[col]]), ]
                       train_data <- train_data[complete.cases(train_data[complete_cols]), ]
                       
                       if (nrow(train_data) <= 1) {
                         warning("No training data available for column: ", col, ". Skipping imputation.")
                         next
                       }
                       else if (nrow(train_data) > 1) {
                         # Build formula
                         formula <- as.formula(paste(col, "~", paste(complete_cols, collapse = "+")))
                         
                         # Train neural network model
                         nn_model <- nnet::nnet(formula, data = train_data, 
                                                size = size, maxit = maxit, trace = FALSE)
                         
                         # Predict missing values
                         missing_indices <- which(is.na(data[[col]]))
                         predict_data <- data[missing_indices, complete_cols, drop = FALSE]
                         
                         # Ensure prediction data has no missing values
                         complete_predict <- complete.cases(predict_data)
                         if (any(complete_predict)) {
                           predictions <- predict(nn_model, 
                                                  predict_data[complete_predict, , drop = FALSE], 
                                                  type = "class")
                           imputed_data[missing_indices[complete_predict], col] <- as.numeric(as.character(predictions))
                         }
                       }
                     }
                     imputed_data
                   },
                   "mi" = {
                     message("Using Multiple Imputation (MICE)...")
                     dots <- list(...)
                     m <- ifelse("m" %in% names(dots), dots$m, 5)
                     maxit <- ifelse("maxit" %in% names(dots), dots$maxit, 10)
                     
                     imp <- mice::mice(data, m = m, maxit = maxit, printFlag = FALSE, ...)
                     mice::complete(imp, 1)  # Return first imputed dataset
                   },
                   "logistic" = {
                     message("Using Logistic Regression imputation...")
                     imputed_data <- data
                     
                     for (col in names(which(colSums(is.na(data)) > 0))) {
                       message("  Imputing column: ", col)
                       
                       # Prepare training data
                       train_data <- data[complete.cases(data[[col]]), ]
                       train_data <- train_data[complete.cases(train_data[complete_cols]), ]
                       
                       if (nrow(train_data) <= 1) {
                         warning("No training data available for column: ", col, ". Skipping imputation.")
                         next
                       }
                       else if (nrow(train_data) > 1) {
                         # Build formula
                         formula <- as.formula(paste(col, "~", paste(complete_cols, collapse = "+")))
                         
                         # Train logistic regression model
                         log_model <- stats::glm(formula, data = train_data, family = binomial)
                         
                         # Predict missing values
                         missing_indices <- which(is.na(data[[col]]))
                         predict_data <- data[missing_indices, complete_cols, drop = FALSE]
                         
                         # Ensure prediction data has no missing values
                         complete_predict <- complete.cases(predict_data)
                         if (any(complete_predict)) {
                           probabilities <- predict(log_model, 
                                                    predict_data[complete_predict, , drop = FALSE], 
                                                    type = "response")
                           predictions <- ifelse(probabilities > 0.5, 1, 0)
                           imputed_data[missing_indices[complete_predict], col] <- predictions
                         }
                       }
                     }
                     imputed_data
                   }
  )
  
  message("Imputation completed using method: ", method)
  message("Original missing values: ", sum(is.na(data)))
  message("Remaining missing values: ", sum(is.na(result)))
  
  return(result)
}
