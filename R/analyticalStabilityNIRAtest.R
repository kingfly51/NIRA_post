#' Bootstrap Stability of the Analytical NIRA Intervention Analysis
#'
#' @description
#' Quantifies the sampling uncertainty of \code{\link{analyticalNIRAtest}} and
#' assesses how stable the ranking of intervention targets is across bootstrap
#' resamples.  Each repetition resamples the original participants with
#' replacement, re-estimates the Ising network, and recomputes the exact
#' analytical solution on the resampled network.  No Monte Carlo samples are
#' simulated, so the only source of variation is sampling error in the
#' estimated network.
#'
#' The output includes the usual bootstrap inference (standard errors,
#' percentile 95 \% confidence intervals, an approximate two-sided
#' \eqn{p}-value for \eqn{H_0: \Delta_i = 0}, and multiplicity-adjusted
#' \eqn{p}-values) together with a \strong{ranking-stability plot}: for every
#' node, the proportion of bootstrap resamples in which it occupied each rank
#' (by absolute intervention effect).  The plot is produced by
#' \code{\link{plotMaxN}} and is the bootstrap analogue of the repeated-
#' simulation stability plot in \code{\link{stabilityNIRAtest}} /
#' \code{\link{findMaxN}}.
#'
#' @param data A binary (0/1) \code{matrix} or \code{data.frame} of dimensions
#'   \eqn{n \times p} (\eqn{n} participants, \eqn{p} nodes).  Rows with missing
#'   values are removed with a warning.  Column names are used as node labels.
#' @param perturbation_type Character string specifying the direction of the
#'   intervention: \code{"aggravating"} (threshold increased, node more likely
#'   to be 1) or \code{"alleviating"} (threshold decreased, node more likely
#'   to be 0).
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
#'   distribution.  Default: \code{2}.
#' @param nReps Positive integer.  Number of bootstrap resamples.
#'   Default: \code{1000}.  Each resample re-estimates the Ising network, so
#'   large values are slow in serial mode; set \code{parallel = TRUE} for
#'   \code{nReps} \eqn{\geq 1000}.
#' @param method Multiple-comparison correction method passed to
#'   \code{\link[stats]{p.adjust}}: \code{"bonferroni"} (default),
#'   \code{"BH"}, \code{"holm"}, \code{"hochberg"}, \code{"hommel"},
#'   \code{"BY"}, \code{"fdr"}, or \code{"none"}.
#' @param seed Integer random seed for reproducibility.  Default: \code{2025}.
#'   Passed to \code{set.seed()} in serial mode and to
#'   \code{parallel::clusterSetRNGStream()} in parallel mode.
#' @param parallel Logical.  If \code{TRUE}, resamples are distributed across
#'   CPU cores.  Default: \code{FALSE}.
#' @param ncores Positive integer or \code{NULL}.  Number of CPU cores for
#'   parallel execution.  If \code{NULL} (default) and
#'   \code{parallel = TRUE}, uses \code{parallel::detectCores() - 1} cores.
#' @param resample_default Character string or \code{NULL} specifying the
#'   estimator used to re-estimate the Ising network on every bootstrap
#'   resample.  The default \code{NULL} matches the point estimate: supplied
#'   parameters imply \code{bootnet::estimateNetwork(default = "IsingFit")};
#'   otherwise the fast saturated pseudo-likelihood estimator is used.  Set
#'   \code{"fast"} to force the fast estimator, or any \code{bootnet}
#'   \code{default} (e.g., \code{"IsingFit"}, \code{"mgm"}) to force a
#'   bootnet estimator.
#' @param nRanks Positive integer or \code{NULL}.  Number of rank positions to
#'   report and plot.  \code{NULL} (default) uses all \eqn{p} positions.
#'   Reduce it for readability in large networks.
#'
#' @section Interpretation:
#' The absolute threshold shift is computed once from the point-estimate
#' thresholds and held fixed across bootstrap resamples.  The intervals thus
#' condition on that intervention magnitude; they do not include uncertainty
#' in the threshold standard deviation used to define it.  Tied ranks share
#' credit equally, so every \code{percenttop_k} column sums to one.
#' Supplied parameters do not identify their estimator or tuning settings;
#' verify that the point fit and bootstrap use compatible estimation settings.
#'
#' @return
#' A named \code{list} with the following elements:
#' \describe{
#'   \item{\code{stat}}{A \code{data.frame} with one row per intervention node
#'     (rows named by node label).  Columns: \code{mean_activation} (exact),
#'     \code{ci_lower}/\code{ci_upper} (percentile 95 \% CI of the resampled
#'     mean activation), \code{effect} (exact), \code{glass_d} (exact),
#'     \code{se_effect}, \code{ci_effect_lower}/\code{ci_effect_upper}
#'     (percentile 95 \% CI of \eqn{\Delta_i}), \code{p},
#'     \code{p.adjust}, \code{stars}, and \code{percenttop_1}…
#'     \code{percenttop_nRanks} (proportion of resamples at each rank).}
#'   \item{\code{exact}}{A \code{data.frame} with the exact mean activation of
#'     the \code{"original"} network and of each one-node intervention.}
#'   \item{\code{resample_effects}}{A numeric \code{matrix} of dimensions
#'     \code{nReps} \eqn{\times} \eqn{p} holding the resampled intervention
#'     effects.}
#'   \item{\code{maxn}}{A \code{data.frame} with \code{percenttop_*}
#'     proportions and \code{repeattop_*} counts for ranks
#'     \eqn{1 \ldots} \code{nRanks}; ready for \code{\link{plotMaxN}}.}
#'   \item{\code{plot}}{A \code{\link[ggplot2]{ggplot}} ranking-stability plot
#'     (stacked bars of rank occupancy across bootstrap resamples).}
#'   \item{\code{parameters}}{A named \code{list} recording the network,
#'     perturbation, resampling and plotting settings used.}
#' }
#'
#' @examples
#' \dontrun{
#' data("single_gds", package = "NIRApost")
#'
#' res <- analyticalStabilityNIRAtest(
#'   data              = single_gds,
#'   perturbation_type = "aggravating",
#'   nReps             = 1000,
#'   parallel          = TRUE,
#'   seed              = 2025
#' )
#' res$stat
#' res$plot            # ranking-stability plot
#' plotMaxN(res$maxn, low = 1, high = 7)
#' }
#'
#' @author
#' Fei Wang \email{bjnwangfei0501@@outlook.com}
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{analyticalNIRAtest}} — exact solution without
#'     confidence intervals.
#'   \item \code{\link{plotMaxN}} / \code{\link{findMaxN}} — the
#'     Monte-Carlo stability plot and ranking table.
#'   \item \code{\link{stabilityNIRAtest}} — repeated-simulation stability
#'     assessment.
#' }
#'
#' @export
analyticalStabilityNIRAtest <- function(data,
                                        perturbation_type = c("aggravating",
                                                              "alleviating"),
                                        edge_weights = NULL,
                                        thresholds = NULL,
                                        amount_of_SDs_perturbation = 2,
                                        nReps = 1000,
                                        method = "bonferroni",
                                        seed = 2025,
                                        parallel = FALSE,
                                        ncores = NULL,
                                        resample_default = NULL,
                                        nRanks = NULL) {

  perturbation_type <- match.arg(perturbation_type)

  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1))
    stop("'seed' must be a single number or NULL.")

  prepared <- analytical_prepare_data(data)
  data       <- prepared$data
  p          <- prepared$p
  node_names <- prepared$node_names

  if (xor(is.null(thresholds), is.null(edge_weights)))
    stop("Supply BOTH 'thresholds' and 'edge_weights', or neither.")

  valid_m <- c("holm", "hochberg", "hommel", "bonferroni", "BH",
               "BY", "fdr", "none")
  if (!method %in% valid_m)
    stop("Invalid method. Choose from: ", paste(valid_m, collapse = ", "))

  if (!is.numeric(nReps) || length(nReps) != 1 || nReps < 1 ||
      nReps != floor(nReps))
    stop("'nReps' must be a single positive integer.")
  nReps <- as.integer(nReps)

  if (!is.null(nRanks)) {
    if (!is.numeric(nRanks) || length(nRanks) != 1 || nRanks < 1 ||
        nRanks != floor(nRanks) || nRanks > p)
      stop("'nRanks' must be NULL or an integer between 1 and p = ", p, ".")
    nRanks <- as.integer(nRanks)
  }
  if (is.null(nRanks)) nRanks <- p

  # Resolve the estimator before fitting: the default depends on whether the
  # user supplied the point-estimate parameters.
  resamp <- analytical_resolve_resample(resample_default, thresholds)

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

  resample_means <- bootstrap_exact_means(
    data          = data,
    delta_spin    = delta_spin,
    source_idx    = idx,
    outcome_idx   = idx,
    nReps         = nReps,
    use_bootnet   = resamp$use_bootnet,
    resample_default = resamp$resample_default,
    seed          = seed,
    parallel      = parallel,
    ncores        = ncores)
  colnames(resample_means) <- c("original", node_names)

  core <- analytical_stability_core(resample_means, effects, ma_nodes,
                                    node_names, sd_original, nRanks, method)

  exact <- data.frame(
    condition       = c("original", node_names),
    mean_activation = c(ma_original, ma_nodes),
    stringsAsFactors = FALSE
  )

  list(
    stat             = core$stat,
    exact            = exact,
    resample_effects = core$resample_effects,
    maxn             = core$maxn,
    plot             = core$plot,
    parameters       = list(
      perturbation_type          = perturbation_type,
      amount_of_SDs_perturbation = amount_of_SDs_perturbation,
      delta                      = delta,
      nReps                      = core$B,
      method                     = method,
      seed                       = seed,
      resample_default           = if (resamp$use_bootnet)
                                     resamp$resample_default else "fast",
      nRanks                     = ncol(core$maxn) / 2L,
      thresholds                 = thresholds,
      edge_weights               = edge_weights
    )
  )
}
