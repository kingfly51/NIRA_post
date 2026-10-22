#' Exact Analytical Bridge (Two-Community) NIRA Intervention Analysis
#'
#' @description
#' Implements the \emph{exact} (analytical) version of a \strong{bridge / dual
#' network} NIRA intervention: for an Ising network partitioned into two
#' communities, it intervenes on each node of a \emph{source} community and
#' computes, \strong{exactly}, the resulting mean number of activated nodes in
#' the \emph{outcome} community.  This is the analytical counterpart of the
#' bridge-network procedure in \code{Bridge_network_NIRA.R}, in which nodes of
#' one community are perturbed and the total score of the other community is
#' the outcome of interest.
#'
#' The exact computation enumerates all \eqn{2^p} activation microstates of
#' the full Ising model (Cui et al., 2026) and tracks the number of active
#' nodes \emph{within the outcome community} (the column subset
#' \eqn{B}), rather than over all \eqn{p} nodes as in
#' \code{\link{analyticalNIRAtest}}.  A one-node intervention changes only the
#' single-node field of the source node, so the perturbed mean is obtained by
#' exactly reweighting the microstates by \eqn{e^{\delta s_i}} — no
#' re-enumeration is needed for the \eqn{p_A} source interventions.  This
#' eliminates Monte Carlo error entirely.
#'
#' This function returns only exact quantities and an exact plot
#' \strong{without confidence intervals}.  Sampling uncertainty and ranking
#' stability are handled by \code{\link{analyticalBridgeStabilityNIRAtest}},
#' which bootstraps the original data, recomputes the exact solution on every
#' resample, and draws a ranking-stability plot analogous to
#' \code{\link{plotMaxN}}.
#'
#' @param data A binary (0/1) \code{matrix} or \code{data.frame} of dimensions
#'   \eqn{n \times p} (\eqn{n} participants, \eqn{p} nodes).  Rows with missing
#'   values are removed with a warning.  Column names are used as node labels.
#' @param groups A character/factor vector of length \eqn{p} giving the
#'   community label of every node (column), e.g. the \code{type} vector
#'   \code{rep(c("JKXW", "MH"), c(8, 3))} used in \code{Bridge_network_NIRA.R}.
#'   Column names of \code{data} are used as node labels within each community.
#'   A single-community partition (all nodes share one label) is a valid
#'   degenerate case: with \code{intervention_group} and \code{outcome_group}
#'   both set to that label it reproduces \code{\link{analyticalNIRAtest}}.
#' @param intervention_group Character string: the label of the \strong{source}
#'   community whose nodes are intervened upon one at a time.
#' @param outcome_group Character string: the label of the \strong{outcome}
#'   community whose mean activation is computed.
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
#'   threshold parameters (e.g., \code{fit$intercepts} from
#'   \code{bootnet::estimateNetwork}).  Must be supplied together with
#'   \code{edge_weights}.  If \code{NULL}, estimated from \code{data}.
#' @param amount_of_SDs_perturbation Positive numeric value specifying the
#'   magnitude of the perturbation in standard deviations of the threshold
#'   distribution.  Default: \code{2} (as in the NIRA tutorial).
#'
#' @section Interpretation:
#' Exact computation does not guarantee a nonzero intervention effect: the
#' reported effects are model-implied differences in outcome-community mean
#' activation, not significance tests.  Conduct sampling-uncertainty and
#' ranking-stability checks with
#' \code{\link{analyticalBridgeStabilityNIRAtest}}.
#'
#' @return
#' A named \code{list} with the following elements:
#' \describe{
#'   \item{\code{stat}}{A \code{data.frame} with one row per \strong{source}
#'     node (rows named by node label) and columns:
#'     \describe{
#'       \item{\code{mean_activation}}{Exact mean number of activated nodes in
#'         the \strong{outcome} community after the intervention.}
#'       \item{\code{effect}}{Intervention effect on the outcome community,
#'         \eqn{\Delta_i = \mathrm{MA}_i - \mathrm{MA}_0} (exact).}
#'       \item{\code{glass_d}}{Standardised effect, i.e. \strong{Glass's
#'         \eqn{\Delta}}: \eqn{\Delta_i / \mathrm{SD}(\mathrm{MA}_0)}, using
#'         the exact standard deviation of the \emph{original} outcome
#'         activation distribution as the standardiser.  This is \emph{not}
#'         Cohen's \eqn{d}, which uses a pooled standard deviation; it
#'         carries the sign of \code{effect} (positive for aggravating,
#'         negative for alleviating).}
#'     }
#'   }
#'   \item{\code{exact}}{A \code{data.frame} with the exact outcome-community
#'     mean activation of the \code{"original"} network and of each one-node
#'     source intervention.}
#'   \item{\code{plot}}{A \code{\link[ggplot2]{ggplot}} object showing the
#'     exact outcome-community mean activation of each condition.  It
#'     deliberately contains no confidence interval or significance
#'     annotation.  Save with \code{ggplot2::ggsave()}.}
#'   \item{\code{parameters}}{A named \code{list} recording the groups,
#'     intervention/outcome communities, thresholds, edge weights and
#'     perturbation definition used.}
#' }
#'
#' @section Method notes:
#' \itemize{
#'   \item The exact computation transforms the \eqn{\{0,1\}}-coded Ising
#'     parameters into the spin (\eqn{\pm 1}) convention and enumerates all
#'     \eqn{2^p} microstates \strong{once per network} with a numerically
#'     stable log-sum-exp scheme.  The outcome statistic is the number of
#'     activated nodes within the outcome community,
#'     \eqn{n_B(s) = (\sum_{j \in B} s_j + p_B) / 2}.  A one-node source
#'     intervention changes only the single-node field, so the perturbed
#'     outcome mean is obtained by exactly reweighting the microstates by
#'     \eqn{e^{\delta s_i}} — no re-enumeration is needed for the \eqn{p_A}
#'     source interventions.
#' }
#'
#' @examples
#' \dontrun{
#' data("single_gds", package = "NIRApost")
#'
#' # Arbitrary two-community partition of the 7 GAD nodes for illustration
#' groups <- c(rep("A", 4), rep("B", 3))
#'
#' res <- analyticalBridgeNIRAtest(
#'   data               = single_gds,
#'   groups             = groups,
#'   intervention_group = "A",     # intervene on community A's nodes
#'   outcome_group      = "B",     # measure community B's mean activation
#'   perturbation_type  = "aggravating"
#' )
#' res$stat            # exact per-source-node effects on the outcome community
#' res$exact           # exact outcome-community mean activation
#' res$plot
#' }
#'
#' @author
#' Fei Wang \email{bjnwangfei0501@@outlook.com}
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{analyticalNIRAtest}} — the single-community (all-node
#'     outcome) analytical workflow; this function is its two-community
#'     generalisation.
#'   \item \code{\link{analyticalBridgeStabilityNIRAtest}} — bootstrap-based
#'     sampling uncertainty and ranking-stability plot.
#'   \item \code{\link{permutationNIRAtest}} — Monte-Carlo-based permutation
#'     test (simulated samples).
#'   \item \code{\link{stabilityNIRAtest}} / \code{\link{findMaxN}} —
#'     Monte-Carlo-based stability assessment.
#'   \item \code{\link[Isinglandr]{make_2d_Isingland}} — the exact
#'     enumeration engine whose output this function reproduces for the
#'     full-network outcome.
#' }
#'
#' @export
analyticalBridgeNIRAtest <- function(data,
                                     groups,
                                     intervention_group,
                                     outcome_group,
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

  # ---- groups / communities ----
  if (length(groups) != p)
    stop("'groups' must have length p = ", p, " (one label per column).")
  groups <- as.character(groups)
  if (anyNA(groups) || any(!nzchar(groups)))
    stop("'groups' must not contain missing or empty labels.")
  # A single-community partition is a valid (degenerate) case: it reduces the
  # analysis to the all-node outcome of analyticalNIRAtest.
  if (!(intervention_group %in% groups))
    stop("'intervention_group' = ", dQuote(intervention_group),
         " is not a label in 'groups'.")
  if (!(outcome_group %in% groups))
    stop("'outcome_group' = ", dQuote(outcome_group),
         " is not a label in 'groups'.")
  source_idx  <- which(groups == intervention_group)
  outcome_idx <- which(groups == outcome_group)
  if (length(source_idx) < 1)
    stop("'intervention_group' contains no nodes.")
  if (length(outcome_idx) < 1)
    stop("'outcome_group' contains no nodes.")
  source_names <- node_names[source_idx]

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

  o <- ising_exact_expectation(thresholds, edge_weights, delta_spin,
                               source_idx = source_idx,
                               outcome_idx = outcome_idx)

  ma_original <- o$mean_original
  sd_original <- o$sd_original
  ma_nodes    <- o$mean_nodes
  effects     <- ma_nodes - ma_original

  stat <- data.frame(
    node            = source_names,
    mean_activation = ma_nodes,
    effect          = effects,
    glass_d         = effects / sd_original,
    row.names       = source_names,
    stringsAsFactors = FALSE
  )

  exact <- data.frame(
    condition       = c("original", source_names),
    mean_activation = c(ma_original, ma_nodes),
    stringsAsFactors = FALSE
  )

  gp <- analytical_exact_plot(
    ma_original, ma_nodes, source_names, perturbation_type,
    ylab = paste0("Mean activation of ", dQuote(outcome_group), " (exact)"))

  list(
    stat       = stat,
    exact      = exact,
    plot       = gp,
    parameters = list(
      perturbation_type            = perturbation_type,
      amount_of_SDs_perturbation   = amount_of_SDs_perturbation,
      delta                        = delta,
      groups                       = groups,
      intervention_group           = intervention_group,
      outcome_group                = outcome_group,
      source_idx                   = source_idx,
      outcome_idx                  = outcome_idx,
      thresholds                   = thresholds,
      edge_weights                 = edge_weights
    )
  )
}
