#' Check Missing Data Patterns and Recommend Imputation Methods
#'
#' @description
#' Diagnoses the missing data structure of a dataset by classifying columns
#' as complete or incomplete, and returns a list of imputation methods
#' appropriate for the detected pattern. This function is the recommended
#' first step in the NIRApost data-preparation pipeline: its output
#' directly guides selection of the \code{method} argument in
#' \code{\link{imputeData}}.
#'
#' When at least one complete column exists, all eight imputation methods
#' are available because feature-based methods (KNN, SVM, decision tree,
#' etc.) can use those complete columns as predictors.  When every column
#' has at least one missing value, only methods that do not require
#' complete predictors are offered (\code{"mode"}, \code{"random_forest"},
#' \code{"mi"}).
#'
#' @param data A \code{data.frame} or \code{matrix} to be examined.
#'   The object must contain named columns. Row and column order is
#'   preserved in the output.
#'
#' @return
#' An invisibly returned named \code{list} with three elements:
#'
#' \describe{
#'   \item{\code{complete_cols}}{Character vector of column names that
#'     contain \strong{no} missing values (\code{NA}).  An empty character
#'     vector when every column has at least one \code{NA}.}
#'   \item{\code{incomplete_cols}}{Character vector of column names that
#'     contain \strong{one or more} missing values.  An empty character
#'     vector when the data are complete.}
#'   \item{\code{available_methods}}{Character vector of imputation method
#'     names that are valid given the missing-data pattern.  Passed
#'     directly to the \code{method} argument of \code{\link{imputeData}}.
#'     Possible values: \code{"mode"}, \code{"knn"}, \code{"random_forest"},
#'     \code{"decision_tree"}, \code{"svm"}, \code{"ann"}, \code{"mi"},
#'     \code{"logistic"}.}
#' }
#'
#' Diagnostic messages are printed to the console; results are
#' returned invisibly so the function may be called for its side-effects
#' only, or its output captured with \code{<-}.
#'
#' @section Method Recommendations:
#' \tabular{ll}{
#'   \strong{Condition} \tab \strong{Available methods} \cr
#'   At least one complete column \tab
#'     \code{mode}, \code{knn}, \code{random_forest}, \code{decision_tree},
#'     \code{svm}, \code{ann}, \code{mi}, \code{logistic} \cr
#'   All columns have missing values \tab
#'     \code{mode}, \code{random_forest}, \code{mi}
#' }
#'
#' \describe{
#'   \item{\code{mode}}{Replaces \code{NA}s with the most frequent
#'     observed value; fastest but ignores relationships between columns.}
#'   \item{\code{knn}}{k-Nearest Neighbours (via \pkg{VIM}); balances
#'     speed and accuracy for moderate-sized datasets.}
#'   \item{\code{random_forest}}{MissForest algorithm (via \pkg{missForest});
#'     handles mixed data types and nonlinear relationships; recommended
#'     when patterns are complex.}
#'   \item{\code{decision_tree}}{CART imputation (via \pkg{rpart}); fast and
#'     interpretable; requires at least one complete predictor column.}
#'   \item{\code{svm}}{Support Vector Machine imputation (via \pkg{e1071});
#'     good for binary outcomes; requires at least one complete column.}
#'   \item{\code{ann}}{Artificial Neural Network imputation (via \pkg{nnet});
#'     captures nonlinear patterns; requires at least one complete column.}
#'   \item{\code{mi}}{Multiple Imputation by Chained Equations
#'     (via \pkg{mice}); the gold-standard approach for binary Ising data;
#'     default method in \code{\link{imputeData}}.}
#'   \item{\code{logistic}}{Logistic regression imputation; appropriate for
#'     binary outcomes; requires at least one complete column.}
#' }
#'
#' @examples
#' ## Example 1: dataset with a mix of complete and incomplete columns
#' df <- data.frame(
#'   A = c(1, 0, 1, 0, 1),         # complete
#'   B = c(NA, 0, NA, 1, 0),       # has NAs
#'   C = c(1, 1, 0, NA, NA)        # has NAs
#' )
#' result <- checkMissing(df)
#' result$complete_cols       # "A"
#' result$incomplete_cols     # c("B", "C")
#' result$available_methods   # all 8 methods
#'
#' ## Example 2: all columns incomplete
#' df_all_na <- data.frame(
#'   X = c(1, NA, 1),
#'   Y = c(NA, 0, 1)
#' )
#' result2 <- checkMissing(df_all_na)
#' result2$available_methods  # "mode", "random_forest", "mi"
#'
#' ## Example 3: use result to drive imputeData()
#' \dontrun{
#' info <- checkMissing(single_gds)
#' if ("mi" %in% info$available_methods) {
#'   clean <- imputeData(single_gds, method = "mi")
#' }
#' }
#'
#' @author
#' Fei Wang \email{bjnwangfei0501@@outlook.com}
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{imputeData}} — performs imputation using the methods
#'     listed by \code{checkMissing}.
#'   \item \code{\link[mice]{mice}} — underlying engine for the \code{"mi"}
#'     method.
#'   \item \code{\link[missForest]{missForest}} — underlying engine for the
#'     \code{"random_forest"} method.
#'   \item \code{\link[VIM]{kNN}} — underlying engine for the \code{"knn"}
#'     method.
#' }
#'
#' @export
checkMissing <- function(data) {

  if (!is.data.frame(data) && !is.matrix(data))
    stop("Input must be a data frame or matrix.")

  complete_cols   <- names(which(colSums(is.na(data)) == 0))
  incomplete_cols <- names(which(colSums(is.na(data)) >  0))

  if (length(complete_cols) > 0) {
    available_methods <- c("mode", "knn", "random_forest", "decision_tree",
                           "svm", "ann", "mi", "logistic")
    message("Complete columns (no NAs): ",
            paste(complete_cols, collapse = ", "))
  } else {
    available_methods <- c("mode", "random_forest", "mi")
    message("All columns contain missing values.")
  }
  message("Available imputation methods: ",
          paste(available_methods, collapse = ", "))

  invisible(list(
    complete_cols     = complete_cols,
    incomplete_cols   = incomplete_cols,
    available_methods = available_methods
  ))
}