#' Impute Missing Values in Binary Ising Network Data
#'
#' @description
#' Provides a unified interface to eight statistical and machine-learning
#' imputation methods for datasets used in Ising network analysis.
#' The function validates method availability for the given data pattern,
#' applies the chosen algorithm, and reports the number of imputed values.
#'
#' Binary (0/1) Ising data require methods that respect the discrete scale.
#' The default method \code{"mi"} (Multiple Imputation by Chained Equations)
#' is the recommended choice for publication-quality analyses because it
#' propagates uncertainty and produces unbiased estimates.  When computational
#' resources are limited, \code{"random_forest"} or \code{"mode"} offer
#' fast, reasonable alternatives.
#'
#' Use \code{\link{checkMissing}} first to inspect the missing-data pattern
#' and confirm which methods are available for your dataset.
#'
#' @param data A \code{data.frame} or \code{matrix} containing the dataset
#'   to be imputed.  Must have named columns.  The object is returned
#'   unchanged if no \code{NA}s are present.
#' @param method Character string selecting the imputation algorithm.
#'   One of \code{"mode"}, \code{"knn"}, \code{"random_forest"},
#'   \code{"decision_tree"}, \code{"svm"}, \code{"ann"}, \code{"mi"}
#'   (default), or \code{"logistic"}.  See the
#'   \strong{Method Details} section for guidance.
#' @param ... Additional arguments forwarded to the underlying imputation
#'   function.  Common options:
#'   \itemize{
#'     \item \code{k} — number of neighbours for \code{"knn"} (default 5).
#'     \item \code{ntree} — number of trees for \code{"random_forest"}
#'       (default 100).
#'     \item \code{maxiter} — maximum iterations for \code{"random_forest"}
#'       (default 10) or \code{"mi"} (default 10).
#'     \item \code{m} — number of imputed datasets for \code{"mi"}
#'       (default 5); the first dataset is returned.
#'     \item \code{size} — hidden units for \code{"ann"} (default 5).
#'     \item \code{maxit} — maximum training epochs for \code{"ann"}
#'       (default 200).
#'   }
#'
#' @return
#' A complete dataset of the same class and dimensions as \code{data},
#' with all \code{NA} values replaced by imputed values.  The original
#' data are returned unchanged if no missing values are detected.
#'
#' The following summary messages are printed:
#' \itemize{
#'   \item Imputation method used.
#'   \item Number of \code{NA}s in the original data.
#'   \item Number of \code{NA}s remaining after imputation (should be 0
#'     for methods that guarantee complete cases; may be non-zero for
#'     column-wise methods when predictor columns are also missing).
#' }
#'
#' @section Method Details:
#' \describe{
#'   \item{\code{"mode"}}{Replaces each \code{NA} with the most frequent
#'     observed value in that column.  Fastest method; ignores
#'     inter-column relationships.  Suitable as a quick baseline or when
#'     missingness is very low (< 5\%).}
#'   \item{\code{"knn"}}{k-Nearest Neighbours imputation via
#'     \code{\link[VIM]{kNN}}.  Borrows information from similar
#'     observations.  Requires at least one complete column.
#'     Default \code{k = 5}.}
#'   \item{\code{"random_forest"}}{MissForest algorithm via
#'     \code{\link[missForest]{missForest}}.  Iteratively trains
#'     random forests on observed data to predict missing values.
#'     Works well for mixed data types and nonlinear patterns.
#'     Does \emph{not} require any complete column.}
#'   \item{\code{"decision_tree"}}{Column-wise CART imputation via
#'     \code{\link[rpart]{rpart}}.  Trains a classification tree on
#'     rows with observed values, then predicts the missing ones.
#'     Requires at least one complete predictor column.}
#'   \item{\code{"svm"}}{Column-wise Support Vector Machine imputation
#'     via \code{\link[e1071]{svm}}.  Suitable for binary outcomes.
#'     Requires at least one complete predictor column.}
#'   \item{\code{"ann"}}{Column-wise Artificial Neural Network
#'     imputation via \code{\link[nnet]{nnet}}.  Captures nonlinear
#'     relationships; slower than tree-based methods.
#'     Requires at least one complete predictor column.}
#'   \item{\code{"mi"}}{Multiple Imputation by Chained Equations via
#'     \code{\link[mice]{mice}}.  Gold-standard approach; propagates
#'     missing-data uncertainty.  Default \code{m = 5} imputed
#'     datasets; the first is returned.  Does \emph{not} require any
#'     complete column.  \strong{Recommended for publication.}}
#'   \item{\code{"logistic"}}{Column-wise logistic regression
#'     imputation via \code{\link[stats]{glm}}.  Appropriate for
#'     binary 0/1 outcomes.  Requires at least one complete predictor
#'     column.}
#' }
#'
#' @section Method Recommendations:
#' \tabular{lll}{
#'   \strong{Scenario} \tab \strong{Recommended method} \tab
#'     \strong{Reason} \cr
#'   Publication analysis \tab \code{"mi"} \tab
#'     Unbiased; uncertainty propagation \cr
#'   Large dataset (> 10 000 rows) \tab \code{"random_forest"} \tab
#'     Scalable; handles mixed types \cr
#'   All columns have NAs \tab \code{"mi"} or \code{"random_forest"} \tab
#'     Do not require complete predictors \cr
#'   Very low missingness (< 5\%) \tab \code{"mode"} \tab
#'     Fast; acceptable bias \cr
#'   Binary outcome, ≥ 1 complete col \tab \code{"logistic"} \tab
#'     Theoretically appropriate for 0/1 \cr
#' }
#'
#' @examples
#' ## Example 1: simple mode imputation
#' df <- data.frame(
#'   A = c(1, 0, 1, 0, 1),
#'   B = c(NA, 0, NA, 1, 0),
#'   C = c(1, 1, 0, NA, NA)
#' )
#' imputed_mode <- imputeData(df, method = "mode")
#' sum(is.na(imputed_mode))  # 0
#'
#' ## Example 2: multiple imputation (recommended)
#' \dontrun{
#' data("single_gds")
#' imputed_mi <- imputeData(single_gds, method = "mi", m = 5, maxit = 10)
#' }
#'
#' ## Example 3: KNN with custom k
#' \dontrun{
#' imputed_knn <- imputeData(df, method = "knn", k = 3)
#' }
#'
#' ## Example 4: random forest with more trees
#' \dontrun{
#' imputed_rf <- imputeData(df, method = "random_forest",
#'                          ntree = 200, maxiter = 15)
#' }
#'
#' ## Example 5: full pipeline with checkMissing
#' \dontrun{
#' info  <- checkMissing(single_gds)
#' clean <- imputeData(single_gds, method = info$available_methods[1])
#' }
#'
#' @author
#' Fei Wang \email{bjnwangfei0501@@outlook.com}
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{checkMissing}} — diagnose missing patterns and select
#'     an appropriate method before calling \code{imputeData}.
#'   \item \code{\link[mice]{mice}} — the underlying engine for
#'     \code{method = "mi"}.
#'   \item \code{\link[missForest]{missForest}} — the underlying engine for
#'     \code{method = "random_forest"}.
#'   \item \code{\link[VIM]{kNN}} — the underlying engine for
#'     \code{method = "knn"}.
#'   \item \code{\link[rpart]{rpart}} — the underlying engine for
#'     \code{method = "decision_tree"}.
#'   \item \code{\link[e1071]{svm}} — the underlying engine for
#'     \code{method = "svm"}.
#'   \item \code{\link[nnet]{nnet}} — the underlying engine for
#'     \code{method = "ann"}.
#' }
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

  if (!is.data.frame(data) && !is.matrix(data))
    stop("Input must be a data frame or matrix.")

  if (sum(is.na(data)) == 0) {
    message("No missing values detected. Returning original data unchanged.")
    return(data)
  }

  valid_methods <- c("mode", "knn", "random_forest", "decision_tree",
                     "svm", "ann", "mi", "logistic")
  if (!method %in% valid_methods)
    stop("Invalid method '", method, "'. Choose from: ",
         paste(valid_methods, collapse = ", "))

  complete_cols <- names(which(colSums(is.na(data)) == 0))
  feature_methods <- c("knn", "decision_tree", "svm", "ann", "logistic")
  if (method %in% feature_methods && length(complete_cols) == 0)
    stop("Method '", method, "' requires at least one complete column. ",
         "Use checkMissing() to verify available methods.")

  result <- switch(method,

    "mode" = {
      message("Applying mode imputation...")
      mode_data <- data
      for (col in names(data)) {
        if (any(is.na(data[[col]]))) {
          mode_val <- names(sort(table(data[[col]]), decreasing = TRUE))[1]
          mode_data[is.na(data[[col]]), col] <- mode_val
        }
      }
      mode_data
    },

    "knn" = {
      message("Applying k-Nearest Neighbours imputation (VIM::kNN)...")
      dots <- list(...)
      k    <- if ("k" %in% names(dots)) dots$k else 5
      res  <- VIM::kNN(data, k = k, ...)
      imp  <- grep("_imp$", names(res), value = TRUE)
      if (length(imp)) res <- res[, !names(res) %in% imp]
      res[, names(data)]
    },

    "random_forest" = {
      message("Applying Random Forest imputation (missForest)...")
      dots    <- list(...)
      maxiter <- if ("maxiter" %in% names(dots)) dots$maxiter else 10
      ntree   <- if ("ntree"   %in% names(dots)) dots$ntree   else 100
      missForest::missForest(data, maxiter = maxiter, ntree = ntree, ...)$ximp
    },

    "decision_tree" = {
      message("Applying Decision Tree imputation (rpart)...")
      out <- data
      for (col in names(which(colSums(is.na(data)) > 0))) {
        message("  Imputing: ", col)
        tr <- data[stats::complete.cases(data[[col]]), ]
        tr <- tr[stats::complete.cases(tr[complete_cols]), ]
        if (nrow(tr) <= 1) { warning("Skip ", col, ": too few training rows."); next }
        fm <- stats::as.formula(paste(col, "~", paste(complete_cols, collapse = "+")))
        m  <- rpart::rpart(fm, data = tr, method = "class")
        ix <- which(is.na(data[[col]]))
        pd <- data[ix, complete_cols, drop = FALSE]
        ok <- stats::complete.cases(pd)
        if (any(ok))
          out[ix[ok], col] <- as.numeric(as.character(
            stats::predict(m, pd[ok, , drop = FALSE], type = "class")))
      }
      out
    },

    "svm" = {
      message("Applying SVM imputation (e1071)...")
      out <- data
      for (col in names(which(colSums(is.na(data)) > 0))) {
        message("  Imputing: ", col)
        tr <- data[stats::complete.cases(data[[col]]), ]
        tr <- tr[stats::complete.cases(tr[complete_cols]), ]
        if (nrow(tr) <= 1) { warning("Skip ", col, ": too few training rows."); next }
        fm <- stats::as.formula(paste(col, "~", paste(complete_cols, collapse = "+")))
        m  <- e1071::svm(fm, data = tr)
        ix <- which(is.na(data[[col]]))
        pd <- data[ix, complete_cols, drop = FALSE]
        ok <- stats::complete.cases(pd)
        if (any(ok))
          out[ix[ok], col] <- as.numeric(as.character(
            stats::predict(m, pd[ok, , drop = FALSE])))
      }
      out
    },

    "ann" = {
      message("Applying Neural Network imputation (nnet)...")
      dots <- list(...)
      sz   <- if ("size"  %in% names(dots)) dots$size  else 5
      mx   <- if ("maxit" %in% names(dots)) dots$maxit else 200
      out  <- data
      for (col in names(which(colSums(is.na(data)) > 0))) {
        message("  Imputing: ", col)
        tr <- data[stats::complete.cases(data[[col]]), ]
        tr <- tr[stats::complete.cases(tr[complete_cols]), ]
        if (nrow(tr) <= 1) { warning("Skip ", col, ": too few training rows."); next }
        fm <- stats::as.formula(paste(col, "~", paste(complete_cols, collapse = "+")))
        m  <- nnet::nnet(fm, data = tr, size = sz, maxit = mx, trace = FALSE)
        ix <- which(is.na(data[[col]]))
        pd <- data[ix, complete_cols, drop = FALSE]
        ok <- stats::complete.cases(pd)
        if (any(ok))
          out[ix[ok], col] <- as.numeric(as.character(
            stats::predict(m, pd[ok, , drop = FALSE], type = "class")))
      }
      out
    },

    "mi" = {
      message("Applying Multiple Imputation by Chained Equations (mice)...")
      dots  <- list(...)
      m     <- if ("m"     %in% names(dots)) dots$m     else 5
      maxit <- if ("maxit" %in% names(dots)) dots$maxit else 10
      imp   <- mice::mice(data, m = m, maxit = maxit, printFlag = FALSE, ...)
      mice::complete(imp, 1)
    },

    "logistic" = {
      message("Applying Logistic Regression imputation...")
      out <- data
      for (col in names(which(colSums(is.na(data)) > 0))) {
        message("  Imputing: ", col)
        tr <- data[stats::complete.cases(data[[col]]), ]
        tr <- tr[stats::complete.cases(tr[complete_cols]), ]
        if (nrow(tr) <= 1) { warning("Skip ", col, ": too few training rows."); next }
        fm <- stats::as.formula(paste(col, "~", paste(complete_cols, collapse = "+")))
        m  <- stats::glm(fm, data = tr, family = stats::binomial)
        ix <- which(is.na(data[[col]]))
        pd <- data[ix, complete_cols, drop = FALSE]
        ok <- stats::complete.cases(pd)
        if (any(ok)) {
          p <- stats::predict(m, pd[ok, , drop = FALSE], type = "response")
          out[ix[ok], col] <- ifelse(p > 0.5, 1, 0)
        }
      }
      out
    }
  )

  message("Imputation complete. Method: ", method)
  message("  Original NAs : ", sum(is.na(data)))
  message("  Remaining NAs: ", sum(is.na(result)))
  result
}