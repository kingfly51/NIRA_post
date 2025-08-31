#' Check Missing Data Patterns and Recommend Imputation Methods
#'
#' @description
#' Diagnoses missing data patterns in a dataset and recommends appropriate
#' imputation methods based on data completeness. Part of the data preprocessing
#' pipeline for statistical analysis and machine learning.
#'
#' @param data A data frame or matrix containing the dataset to be checked.
#'
#' @return
#' A list with three components:
#'
#' **1. `complete_cols`**
#' Character vector of column names with no missing values
#'
#' **2. `incomplete_cols`**
#' Character vector of column names containing missing values
#'
#' **3. `available_methods`**
#' Character vector of recommended imputation methods suitable for the data pattern
#'
#' @examples
#' ```r
#' # Create dataset with missing patterns
#' df <- data.frame(
#'   A = c(1, 2, 3, 4, 5),
#'   B = c(NA, 2, NA, 4, 5),
#'   C = c(1, 2, 3, NA, NA)
#' )
#'
#' # Check missing patterns
#' result <- check_missing(df)
#'
#' # Access results
#' complete_cols <- result$complete_cols
#' available_methods <- result$available_methods
#' ```
#'
#' @section Method Recommendations:
#' - **Complete columns**: All methods available
#' - **All columns incomplete**: Restricted to robust methods (`mode`, `random_forest`, `mi`)
#'
#' @author
#' Fei Wang,State Key Laboratory of Cognitive Science and Mental Health, Institute of Psychology, Chinese Academy of Sciences (2025-08-31)
#'
#' @seealso
#' * [`mice::mice()`] for multiple imputation implementation
#' * [`missForest::missForest()`] for random forest imputation
#' * [`VIM::kNN()`] for k-nearest neighbors imputation
#'
#' @export
checkMissing <- function(data) {
  # Input validation
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("Input must be a data frame or matrix")
  }

  # Calculate missing patterns
  complete_cols <- names(which(colSums(is.na(data)) == 0))
  incomplete_cols <- names(which(colSums(is.na(data)) > 0))

  # Determine available methods based on missing pattern
  if (length(complete_cols) > 0) {
    available_methods <- c("mode", "knn", "random_forest", "decision_tree",
                           "svm", "ann", "mi", "logistic")
    message("The dataset contains columns with no missing values: ",
            paste(complete_cols, collapse = ", "))
    message("The following methods can be used for imputation: ",
            paste(available_methods, collapse = ", "))
  } else {
    available_methods <- c("mode", "random_forest", "mi")
    message("All columns in the dataset contain missing values")
    message("The following methods can be used for imputation: ",
            paste(available_methods, collapse = ", "))
  }

  # Return structured results
  return(list(
    complete_cols = complete_cols,
    incomplete_cols = incomplete_cols,
    available_methods = available_methods
  ))
}

