#' Batch Permutation Significance Tests Across All NIRA Nodes
#'
#' @description
#' Applies \code{\link{permutationTest}} to every intervention node in a
#' NIRA result and corrects for multiple comparisons, producing both a
#' statistics table and a plotting-ready data frame.
#'
#' Visual inspection of \code{nodeIdentifyR::plotSumScores} alone has four
#' critical limitations that this function addresses:
#' \enumerate{
#'   \item \strong{No formal evidence} — visual comparison cannot confirm
#'     statistical significance.
#'   \item \strong{Unreliable CI overlap} — assessing significance by
#'     checking whether 95 \% error bars overlap is known to be inaccurate,
#'     especially when means are close.
#'   \item \strong{Magnitude sensitivity} — a 2-SD perturbation typically
#'     produces large, visually obvious differences; smaller perturbations
#'     make visual assessment unreliable.
#'   \item \strong{Multiple-comparison inflation} — \eqn{p} simultaneous
#'     comparisons are made with no correction, inflating the family-wise
#'     error rate.
#' }
#'
#' The function internally calls \code{\link{permutationTest}} for each
#' node (5 000 permutations each), assembles the raw p-values, applies the
#' chosen correction, and returns both per-node statistics and a structured
#' data frame ready for \code{\link{plotMeanNIRA}}.
#'
#' @param gs_sumIsingSamplesLong A long-format \code{data.frame} produced
#'   by \code{nodeIdentifyR::prepareDFforPlottingAndANOVA()}.  Must contain
#'   at least two columns:
#'   \describe{
#'     \item{\code{sample}}{Character/factor column identifying each network
#'       condition: \code{"original"} for the unperturbed network,
#'       node names for each intervention condition.}
#'     \item{\code{sumscore}}{Numeric column of simulated total scores.}
#'   }
#' @param method Multiple-comparison correction method passed to
#'   \code{\link[stats]{p.adjust}}.  One of:
#'   \code{"bonferroni"} (default), \code{"BH"}, \code{"holm"},
#'   \code{"hochberg"}, \code{"hommel"}, \code{"BY"}, \code{"fdr"},
#'   \code{"none"}.
#'
#' @return
#' A named \code{list} with two elements:
#'
#' \describe{
#'   \item{\code{stat}}{A \code{data.frame} with one row per intervention
#'     node (rows named by node label) and nine columns:
#'     \describe{
#'       \item{\code{mean_other}}{Mean sum score of the intervention sample.}
#'       \item{\code{sd_other}}{Standard deviation.}
#'       \item{\code{se_other}}{Standard error.}
#'       \item{\code{ci_other}}{95 \% CI half-width (\eqn{1.96 \times \mathrm{SE}}).}
#'       \item{\code{ciLower_other}}{Lower 95 \% CI bound.}
#'       \item{\code{ciUpper_other}}{Upper 95 \% CI bound.}
#'       \item{\code{cohen_d}}{Cohen's \eqn{d} versus the original sample.}
#'       \item{\code{p}}{Raw two-sided permutation p-value.}
#'       \item{\code{p.adjust}}{Adjusted p-value after multiple-comparison
#'         correction.}
#'     }
#'   }
#'   \item{\code{plot_data}}{A \code{data.frame} with one row per network
#'     condition (including the original) and four columns:
#'     \describe{
#'       \item{\code{mean}}{Mean sum score.}
#'       \item{\code{ciLower}}{Lower 95 \% CI bound.}
#'       \item{\code{ciUpper}}{Upper 95 \% CI bound.}
#'       \item{\code{id}}{Integer position index (0 = original).}
#'       \item{\code{node}}{Condition name.}
#'     }
#'     This data frame is passed directly to \code{\link{plotMeanNIRA}}.
#'   }
#' }
#'
#' @section Method Recommendations:
#' \tabular{lll}{
#'   \strong{Scenario} \tab \strong{method} \tab \strong{Rationale} \cr
#'   Standard (recommended) \tab \code{"bonferroni"} \tab
#'     Conservative; controls FWER \cr
#'   Exploratory / many nodes \tab \code{"BH"} or \code{"fdr"} \tab
#'     Controls FDR; more power \cr
#'   Small network (≤ 5 nodes) \tab \code{"holm"} \tab
#'     Uniformly more powerful than Bonferroni \cr
#'   No correction needed \tab \code{"none"} \tab
#'     When only one node is of interest \cr
#' }
#'
#' @examples
#' ## Example 1: standard workflow
#' \dontrun{
#' data("single_gds")
#' gs_fit     <- bootnet::estimateNetwork(single_gds,
#'                 default = "IsingFit", rule = "AND")
#' gs_samples <- nodeIdentifyR::simulateResponses(
#'                 gs_fit$graph, gs_fit$intercepts, "alleviating", 2)
#' gs_long    <- nodeIdentifyR::prepareDFforPlottingAndANOVA(
#'                 nodeIdentifyR::calculateSumScores(gs_samples))
#'
#' stat_result <- permutationNIRAtest(gs_long, method = "bonferroni")
#' head(stat_result$stat)         # per-node statistics
#' plotMeanNIRA(stat_result, "alleviating")
#' }
#'
#' ## Example 2: FDR correction for a large network
#' \dontrun{
#' stat_fdr <- permutationNIRAtest(gs_long, method = "BH")
#' # Nodes significant after FDR correction
#' subset(stat_fdr$stat, p.adjust < 0.05)
#' }
#'
#' ## Example 3: extract effect sizes
#' \dontrun{
#' stat_result <- permutationNIRAtest(gs_long)
#' stat_result$stat[order(abs(stat_result$stat$cohen_d),
#'                        decreasing = TRUE), "cohen_d", drop = FALSE]
#' }
#'
#' @author
#' Fei Wang \email{bjnwangfei0501@@outlook.com}
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{permutationTest}} — the single-comparison engine
#'     called internally for each node.
#'   \item \code{\link{plotMeanNIRA}} — pass \code{stat_result} here to
#'     generate the significance plot.
#'   \item \code{\link[stats]{p.adjust}} — the correction function applied
#'     to the raw p-values.
#'   \item \code{\link{runMgmmAnalysis}} — run this first to confirm NIRA's
#'     theoretical assumptions hold.
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by mutate select across row_number
#' @importFrom tidyr pivot_wider
#' @importFrom stats p.adjust
#' @importFrom rlang .data
#' @export
permutationNIRAtest <- function(gs_sumIsingSamplesLong, method = "bonferroni") {

  if (!is.data.frame(gs_sumIsingSamplesLong))
    stop("Input must be a data.frame from ",
         "nodeIdentifyR::prepareDFforPlottingAndANOVA().")
  if (!all(c("sample", "sumscore") %in% colnames(gs_sumIsingSamplesLong)))
    stop("Input must contain columns 'sample' and 'sumscore'.")
  valid_m <- c("holm", "hochberg", "hommel", "bonferroni",
               "BH", "BY", "fdr", "none")
  if (!method %in% valid_m)
    stop("Invalid method. Choose from: ", paste(valid_m, collapse = ", "))

  # Pivot to wide format
  wide <- gs_sumIsingSamplesLong %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(row = dplyr::row_number()) %>%
    tidyr::pivot_wider(names_from = sample, values_from = sumscore) %>%
    dplyr::select(-row) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric))

  cats  <- unique(gs_sumIsingSamplesLong[[2]])
  orig  <- as.matrix(wide[, "original"])
  other <- as.matrix(wide[, !cats %in% "original"])

  set.seed(2025)
  sl <- as.data.frame(NULL)
  for (i in seq_len(ncol(other))) {
    pl          <- permutationTest(orig, other[, i], 5000)
    sl[i, 1]   <- pl$mean_other
    sl[i, 2]   <- pl$sd_other
    sl[i, 3]   <- pl$se_other
    sl[i, 4]   <- pl$ci_other
    sl[i, 5]   <- pl$ciLower_other
    sl[i, 6]   <- pl$ciUpper_other
    sl[i, 7]   <- pl$cohens_d
    sl[i, 8]   <- pl$p_value
  }
  rownames(sl) <- colnames(other)
  colnames(sl) <- c("mean_other", "sd_other", "se_other", "ci_other",
                    "ciLower_other", "ciUpper_other", "cohen_d", "p")
  sl$p.adjust  <- stats::p.adjust(sl$p, method = method)

  # Build plot_data
  d1           <- sl[, c(1, 5, 6)]
  colnames(d1) <- c("mean", "ciLower", "ciUpper")
  d2           <- cbind(pl$mean_original,
                        pl$ciLower_original,
                        pl$ciUpper_original)
  colnames(d2) <- c("mean", "ciLower", "ciUpper")
  rownames(d2) <- "original"
  d3           <- rbind(d2, d1)
  d3$id        <- 0:(nrow(d3) - 1)
  d3$node      <- rownames(d3)

  list(stat = sl, plot_data = d3)
}