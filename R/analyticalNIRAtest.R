#' Exact Analytical NIRA Intervention Analysis
#'
#' @description
#' Implements the \emph{exact} (analytical) version of the NIRA intervention
#' analysis.  Instead of generating Monte Carlo samples with
#' \code{nodeIdentifyR::simulateResponses}, the mean number of activated nodes
#' is computed \strong{exactly} by enumerating all \eqn{2^p} activation
#' microstates of the Ising model (Cui et al., 2026).  The analytical solution
#' eliminates Monte Carlo error, so the returned quantities are deterministic
#' functions of the fitted network parameters.
#'
#' This function returns only exact quantities and an exact plot
#' \strong{without confidence intervals}.  Uncertainty from estimating the
#' network on a finite sample is deliberately kept in a separate function,
#' \code{\link{analyticalStabilityNIRAtest}}, which bootstraps the original
#' data, recomputes the exact solution on every resample, and draws a
#' ranking-stability plot analogous to \code{\link{plotMaxN}}.
#'
#' @param data A binary (0/1) \code{matrix} or \code{data.frame} of dimensions
#'   \eqn{n \times p} (\eqn{n} participants, \eqn{p} nodes).  Rows with missing
#'   values are removed with a warning.  Column names are used as node labels.
#' @param perturbation_type Character string specifying the direction of the
#'   intervention, matching \code{nodeIdentifyR::simulateResponses}:
#'   \code{"aggravating"} (threshold increased, node more likely to be 1) or
#'   \code{"alleviating"} (threshold decreased, node more likely to be 0).
#' @param edge_weights Optional symmetric \eqn{p \times p} numeric matrix of
#'   Ising edge weights (e.g., \code{fit$graph} from
#'   \code{bootnet::estimateNetwork}).  Must be supplied together with
#'   \code{thresholds}.  If \code{NULL}, both parameters are estimated from
#'   \code{data} with a fast saturated pseudo-likelihood estimator.
#' @param thresholds Optional numeric vector of length \eqn{p} of Ising
#'   threshold parameters (e.g., \code{fit$results$thresholds} from
#'   \code{bootnet::estimateNetwork}).  Must be supplied together with
#'   \code{edge_weights}.  If \code{NULL}, estimated from \code{data}.
#' @param amount_of_SDs_perturbation Positive numeric value specifying the
#'   magnitude of the perturbation in standard deviations of the threshold
#'   distribution.  Default: \code{2} (as in the NIRA tutorial).
#'
#' @section Interpretation:
#' Exact computation does not guarantee a nonzero intervention effect: the
#' reported effects are model-implied differences in mean activation, not
#' significance tests.  Conduct sampling-uncertainty and ranking-stability
#' checks with \code{\link{analyticalStabilityNIRAtest}}.
#'
#' @return
#' A named \code{list} with the following elements:
#' \describe{
#'   \item{\code{stat}}{A \code{data.frame} with one row per intervention node
#'     (rows named by node label) and columns:
#'     \describe{
#'       \item{\code{mean_activation}}{Exact mean number of activated nodes
#'         after the intervention.}
#'       \item{\code{effect}}{Intervention effect
#'         \eqn{\Delta_i = \mathrm{MA}_i - \mathrm{MA}_0} (exact).}
#'       \item{\code{glass_d}}{Standardised effect, i.e. \strong{Glass's
#'         \eqn{\Delta}}: \eqn{\Delta_i / \mathrm{SD}(\mathrm{MA}_0)}, using
#'         the exact standard deviation of the \emph{original} activation
#'         distribution as the standardiser.  This is \emph{not} Cohen's
#'         \eqn{d}, which uses a pooled standard deviation; it carries the
#'         sign of \code{effect} (positive for aggravating, negative for
#'         alleviating).}
#'     }
#'   }
#'   \item{\code{exact}}{A \code{data.frame} with the exact mean activation of
#'     the \code{"original"} network and of each one-node intervention.}
#'   \item{\code{plot}}{A \code{\link[ggplot2]{ggplot}} object showing the
#'     exact mean activation of each condition.  It deliberately contains no
#'     confidence interval or significance annotation.  Save with
#'     \code{ggplot2::ggsave()}.}
#'   \item{\code{parameters}}{A named \code{list} recording the thresholds,
#'     edge weights and perturbation definition used.}
#' }
#'
#' @section Method notes:
#' \itemize{
#'   \item The exact computation transforms the \eqn{\{0,1\}}-coded Ising
#'     parameters into the spin (\eqn{\pm 1}) convention and enumerates all
#'     \eqn{2^p} microstates \strong{once per network} with a numerically
#'     stable log-sum-exp scheme.  A one-node intervention changes only the
#'     single-node field, so the perturbed mean activation is obtained by
#'     exactly reweighting the microstates by \eqn{e^{\Delta_i v_i}} — no
#'     re-enumeration of the \eqn{2^p} states is needed for the \eqn{p}
#'     interventions.  The result is identical to
#'     \code{Isinglandr::make_2d_Isingland(transform = TRUE)}.
#' }
#'
#' @examples
#' ## Exact analysis from the built-in data (fast estimator)
#' \dontrun{
#' data("single_gds", package = "NIRApost")
#'
#' res <- analyticalNIRAtest(
#'   data              = single_gds,
#'   perturbation_type = "aggravating"
#' )
#' res$stat
#' res$plot
#' }
#'
#' ## Reproduce the paper's exact Fig. 5 point estimates by supplying the
#' ## bootnet IsingFit parameters
#' \dontrun{
#' fit <- bootnet::estimateNetwork(single_gds, default = "IsingFit",
#'                                 rule = "AND")
#' res <- analyticalNIRAtest(
#'   data              = single_gds,
#'   edge_weights      = fit$results$weiadj,
#'   thresholds        = fit$results$thresholds,
#'   perturbation_type = "aggravating"
#' )
#' res$exact
#' }
#'
#' @author
#' Fei Wang \email{bjnwangfei0501@@outlook.com}
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{analyticalStabilityNIRAtest}} — bootstrap-based
#'     sampling uncertainty and ranking-stability plot.
#'   \item \code{\link{permutationNIRAtest}} — Monte-Carlo-based permutation
#'     test (simulated samples).
#'   \item \code{\link{stabilityNIRAtest}} / \code{\link{findMaxN}} —
#'     Monte-Carlo-based stability assessment.
#'   \item \code{\link[Isinglandr]{make_2d_Isingland}} — the exact
#'     enumeration engine whose output this function reproduces.
#' }
#'
#' @export
analyticalNIRAtest <- function(data,
                               perturbation_type = c("aggravating",
                                                     "alleviating"),
                               edge_weights = NULL,
                               thresholds = NULL,
                               amount_of_SDs_perturbation = 2) {

  perturbation_type <- match.arg(perturbation_type)

  prepared <- analytical_prepare_data(data)
  data       <- prepared$data
  p          <- prepared$p
  node_names <- prepared$node_names

  if (xor(is.null(thresholds), is.null(edge_weights)))
    stop("Supply BOTH 'thresholds' and 'edge_weights', or neither.")

  net <- analytical_fit_network(data, thresholds, edge_weights, p)
  thresholds   <- net$thresholds
  edge_weights <- net$edge_weights

  dpar       <- analytical_resolve_delta(thresholds,
                                         amount_of_SDs_perturbation,
                                         perturbation_type)
  delta      <- dpar$delta
  delta_spin <- dpar$delta_spin

  idx <- seq_len(p)
  o   <- ising_exact_expectation(thresholds, edge_weights, delta_spin,
                                 source_idx = idx, outcome_idx = idx)

  ma_original <- o$mean_original
  sd_original <- o$sd_original
  ma_nodes    <- o$mean_nodes
  effects     <- ma_nodes - ma_original

  stat <- data.frame(
    node            = node_names,
    mean_activation = ma_nodes,
    effect          = effects,
    glass_d         = effects / sd_original,
    row.names       = node_names,
    stringsAsFactors = FALSE
  )

  exact <- data.frame(
    condition       = c("original", node_names),
    mean_activation = c(ma_original, ma_nodes),
    stringsAsFactors = FALSE
  )

  gp <- analytical_exact_plot(ma_original, ma_nodes, node_names,
                              perturbation_type,
                              ylab = "Mean activation (exact)")

  list(
    stat       = stat,
    exact      = exact,
    plot       = gp,
    parameters = list(
      perturbation_type          = perturbation_type,
      amount_of_SDs_perturbation = amount_of_SDs_perturbation,
      delta                      = delta,
      thresholds                 = thresholds,
      edge_weights               = edge_weights
    )
  )
}
