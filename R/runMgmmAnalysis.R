#' Moderated Graphical Model Analysis — NIRA Prerequisite Test
#'
#' @description
#' Tests whether any node in an Ising network acts as a \emph{moderator}
#' of other pairwise edges.  This test \strong{must be run before} the NIRA
#' simulation because NIRA assumes that perturbing one node's threshold
#' leaves all edge weights unchanged.  If node \eqn{A} moderates the
#' \eqn{B}–\eqn{C} relationship, then altering \eqn{A}'s intercept
#' implicitly changes the \eqn{B}–\eqn{C} edge weight, making it
#' theoretically incorrect to reuse the original edge weight matrix in
#' post-intervention simulations.
#'
#' The function iteratively treats each of the \eqn{p} nodes as a
#' candidate moderator.  For each candidate it:
#' \enumerate{
#'   \item Fits a Moderated Mixed Graphical Model (MGMM) via
#'     \code{\link[mgm]{mgm}}, testing whether that node moderates every
#'     pairwise edge among the remaining nodes.
#'   \item Applies case-resampling via \code{\link[mgm]{resample}} to
#'     obtain resampling-based confidence intervals for each moderation weight.
#'   \item Retains only effects whose 95 \% CI excludes zero (i.e., both
#'     quartile bounds have the same sign).
#' }
#'
#' An \strong{empty} \code{significant_moderators} list confirms that the
#' NIRA assumption holds and the simulation can proceed.  A
#' \strong{non-empty} list identifies which nodes act as moderators and
#' which edges they affect; these should be reported and the NIRA results
#' interpreted with caution.
#'
#' @param data A numeric \code{matrix} or \code{data.frame} of dimensions
#'   \eqn{n \times p} where \eqn{n} is the number of participants and
#'   \eqn{p} is the number of binary (0/1) nodes.
#' @param plotResults Logical.  If \code{TRUE}, a PNG diagnostic plot is
#'   saved to the current working directory for each moderator node
#'   (filename: \code{Moderator_<i>_Plot.png}).  Default: \code{FALSE}.
#' @param rule Character string specifying the edge-selection rule:
#'   \code{"AND"} (default, more conservative — both regression directions
#'   must be significant) or \code{"OR"} (more inclusive — either direction
#'   sufficient).  Consistent with the rule used in
#'   \code{bootnet::estimateNetwork}.
#' @param lambdaGam Numeric EBIC tuning parameter in \eqn{[0, 1]}.
#'   Larger values impose stronger sparsity.  Default: \code{0.25}.
#'   \itemize{
#'     \item \code{0}: BIC-like penalty (less sparse).
#'     \item \code{0.5}: recommended for exploratory analyses.
#'     \item \code{1}: maximum EBIC penalty (most sparse).
#'   }
#' @param nB Integer number of resampling iterations.  Default: \code{100}.
#'   \itemize{
#'     \item \code{nB >= 100}: acceptable for exploratory analyses.
#'     \item \code{nB >= 1000}: recommended for publication.
#'   }
#'   Very small values (\code{nB < 10}) trigger a warning.
#'
#' @return
#' A named \code{list} with two elements:
#'
#' \describe{
#'   \item{\code{all_results}}{A \code{list} of length \eqn{p} (one element
#'     per node tested as moderator).  Each element is itself a
#'     \code{list} containing:
#'     \describe{
#'       \item{\code{model}}{The fitted \code{mgm} object.}
#'       \item{\code{interactions}}{Interaction indicator matrix from
#'         \code{model$interactions$indicator[[2]]}.}
#'       \item{\code{resamplingResult}}{The resampling result object from
#'         \code{\link[mgm]{resample}}.}
#'       \item{\code{table}}{Summary table of edge weights and 95 \% CI
#'         bounds (columns 7–8 removed for brevity).}
#'     }
#'   }
#'   \item{\code{significant_moderators}}{A \code{list} containing only
#'     those nodes for which at least one edge's 95 \% CI excludes zero.
#'     \strong{An empty list means no significant moderation was detected
#'     and NIRA can proceed.}  Each element contains:
#'     \describe{
#'       \item{\code{node}}{Integer index of the moderating node.}
#'       \item{\code{significant_edges}}{Row indices of significantly
#'         moderated edges in \code{table}.}
#'       \item{\code{interaction_details}}{Sub-table of the significant
#'         rows, including edge identifiers, mean weights, and CI bounds.}
#'     }
#'   }
#' }
#'
#' @section Method Recommendations:
#' \tabular{lll}{
#'   \strong{Goal} \tab \strong{Setting} \tab \strong{Rationale} \cr
#'   Fast check \tab \code{nB = 100}, \code{rule = "AND"} \tab
#'     Reasonable for exploratory work \cr
#'   Publication \tab \code{nB = 1000}, \code{rule = "AND"} \tab
#'     Stable CI estimates \cr
#'   Liberal screening \tab \code{rule = "OR"} \tab
#'     Catches more potential moderators \cr
#'   Dense network \tab higher \code{lambdaGam} \tab
#'     More aggressive sparsity control \cr
#' }
#'
#' @examples
#' ## Example 1: quick check on the built-in GAD-7 dataset
#' \dontrun{
#' data("single_gds")
#' sig_moder <- runMgmmAnalysis(
#'   data         = as.matrix(single_gds),
#'   plotResults  = FALSE,
#'   rule         = "AND",
#'   lambdaGam    = 0.25,
#'   nB           = 10        # increase to >= 100 for real analyses
#' )
#' # Check outcome: empty list => NIRA is appropriate
#' print(sig_moder$significant_moderators)
#' }
#'
#' ## Example 2: inspect a specific moderator's resampling table
#' \dontrun{
#' sig_moder <- runMgmmAnalysis(as.matrix(single_gds), nB = 100)
#' # First moderator's resampling table
#' sig_moder$all_results[["Moderator_1"]]$table
#' }
#'
#' ## Example 3: save diagnostic plots and use OR rule
#' \dontrun{
#' sig_moder <- runMgmmAnalysis(
#'   data        = as.matrix(PHQ9),
#'   plotResults = TRUE,      # saves Moderator_*_Plot.png
#'   rule        = "OR",
#'   nB          = 500
#' )
#' }
#'
#' @author
#' Fei Wang \email{bjnwangfei0501@@outlook.com}
#'
#' @seealso
#' \itemize{
#'   \item \code{\link[mgm]{mgm}} — fits the Moderated Mixed Graphical
#'     Model underlying this function.
#'   \item \code{\link[mgm]{resample}} — case-resampling routine used to
#'     compute resampling-based CIs.
#'   \item \code{\link{stabilityNIRAtest}} — run this \emph{after} confirming
#'     no significant moderation effects.
#'   \item \code{\link{permutationNIRAtest}} — significance testing of
#'     NIRA intervention effects.
#' }
#'
#' @importFrom mgm mgm resample plotRes
#' @importFrom grDevices png dev.off
#' @importFrom utils capture.output
#' @export
runMgmmAnalysis <- function(data, plotResults = FALSE, rule = "AND",
                            lambdaGam = 0.25, nB = 100) {

  if (!is.matrix(data) && !is.data.frame(data))
    stop("data must be a matrix or data frame.")
  if (!rule %in% c("AND", "OR"))
    stop("rule must be 'AND' or 'OR'.")
  if (lambdaGam < 0 || lambdaGam > 1)
    stop("lambdaGam must be in [0, 1].")
  if (nB < 10)
    warning("nB < 10 may yield unstable resampling-based confidence intervals.")

  p    <- ncol(data)
  rl   <- vector("list", p)
  names(rl) <- paste0("Moderator_", seq_len(p))

  for (i in seq_len(p)) {
    # Set a deterministic per-node seed before mgm::resample() to ensure
    # reproducibility. mgm::resample() internally calls set.seed(1),
    # which would otherwise override any externally set seed.
    # Using node index i as offset ensures independence across nodes.
    set.seed(1314L + i)
    cat("\n------ Testing node", i, "as moderator ------\n")

    mm <- mgm::mgm(
      data       = as.matrix(data),
      type       = rep("c", p),
      level      = rep(2L, p),
      lambdaSel  = "EBIC",
      lambdaGam  = lambdaGam,
      ruleReg    = rule,
      moderators = i,
      scale      = FALSE
    )

    br  <- mgm::resample(object = mm, data = as.matrix(data), nB = nB)
    tab <- mgm::plotRes(object = br, table = TRUE)[, -c(7, 8)]

    rl[[i]] <- list(
      model        = mm,
      interactions = mm$interactions$indicator[[2]],
      resamplingResult = br,
      table        = tab
    )

    if (plotResults) {
      grDevices::png(paste0("Moderator_", i, "_Plot.png"),
                    width = 800, height = 600)
      mgm::plotRes(object = br)
      grDevices::dev.off()
    }
  }

  # Identify significant moderators (95 % CI excludes zero)
  sm <- list()
  for (i in seq_len(p)) {
    tb  <- rl[[i]]$table
    sig <- which(tb[, "Mod_qtl_low"] > 0 | tb[, "Mod_qtl_high"] < 0)
    if (length(sig) > 0)
      sm[[paste0("Moderator_", i)]] <- list(
        node                = i,
        significant_edges   = sig,
        interaction_details = rl[[i]]$table[sig, ]
      )
  }

  list(all_results = rl, significant_moderators = sm)
}
