#' Repeated-Simulation Stability Test for NIRA
#'
#' @description
#' Assesses whether NIRA results are stable across different random seeds by
#' independently repeating the full simulation procedure \code{nReps} times.
#' Each repetition uses the original Ising network parameters (edge weight
#' matrix and threshold vector) to generate a completely fresh set of
#' simulated samples.  The observed data are \strong{never resampled};
#' this is \emph{repeated simulation} (re-running the simulation from
#' scratch using fixed network parameters), not resampling of the
#' original data.
#'
#' \strong{Why stability testing is necessary:}
#' A single simulation run may give different results for different random
#' seeds — a different node may appear as the "best" intervention target.
#' Only results that remain consistent across many independent runs can be
#' trusted.  If conclusions fluctuate unpredictably across seeds, the NIRA
#' findings cannot be considered reliable.
#'
#' \strong{Note:} Each repetition independently re-runs
#' \code{nodeIdentifyR::simulateResponses} from scratch using the fixed
#' network parameters (edge weights and thresholds) — the observed data
#' are never touched.
#'
#' @param edge_weights Symmetric numeric matrix of dimensions
#'   \eqn{p \times p} containing Ising edge weights.  Typically
#'   \code{fit$graph} from \code{bootnet::estimateNetwork}.  The diagonal
#'   should be zero.
#' @param thresholds Numeric vector of length \eqn{p} containing node
#'   threshold parameters (intercepts from the Ising logistic regression).
#'   Typically \code{fit$intercepts}.
#' @param perturbation_type Character string specifying the direction of
#'   intervention: \code{"alleviating"} (threshold decreased, node more
#'   likely to be 0) or \code{"aggravating"} (threshold increased, node
#'   more likely to be 1).  Must match the type used in
#'   \code{nodeIdentifyR::simulateResponses}.
#' @param amount_of_SDs_perturbation Positive numeric value specifying the
#'   magnitude of the perturbation in standard deviations of the threshold
#'   distribution.  Typical range: 1–3.  Default value used in the NIRA
#'   tutorial: \code{2}.
#' @param nReps Positive integer.  Number of independent simulation
#'   repetitions.  Default: \code{1000}.
#'   \itemize{
#'     \item \code{nReps >= 100}: acceptable for exploratory work.
#'     \item \code{nReps >= 1000}: recommended for publication.
#'   }
#' @param parallel Logical.  If \code{TRUE}, uses \pkg{parallel} to
#'   distribute repetitions across CPU cores.  Default: \code{FALSE}.
#'   Recommended when \eqn{p > 10} or \code{nReps > 100}.
#' @param ncores Positive integer or \code{NULL}.  Number of CPU cores for
#'   parallel execution.  If \code{NULL} (default) and
#'   \code{parallel = TRUE}, the function uses
#'   \code{parallel::detectCores() - 1} cores.  Ignored when
#'   \code{parallel = FALSE}.
#' @param seed Integer random seed for reproducibility.  Default:
#'   \code{2025}.  Passed to \code{set.seed()} in serial mode and to
#'   \code{parallel::clusterSetRNGStream()} in parallel mode.
#'
#' @return
#' A named \code{list} with three elements:
#'
#' \describe{
#'   \item{\code{SimSamples}}{A \code{list} of length \code{nReps}.
#'     Each element is the long-format \code{data.frame} produced by
#'     \code{nodeIdentifyR::prepareDFforPlottingAndANOVA} for one
#'     simulation repetition, with columns \code{sample} (condition name)
#'     and \code{sumscore} (simulated total score).}
#'   \item{\code{mean}}{Numeric \code{matrix} of dimensions
#'     \code{nReps} \eqn{\times} (N nodes + 1).  Columns correspond to
#'     the original condition and each intervention node; rows correspond
#'     to repetitions.  Column names are condition labels.}
#'   \item{\code{sd}}{Numeric \code{matrix} of the same dimensions as
#'     \code{mean}, containing per-repetition within-condition standard
#'     deviations.}
#' }
#'
#' @section Method Recommendations:
#' \tabular{lll}{
#'   \strong{Scenario} \tab \strong{Setting} \tab \strong{Note} \cr
#'   Quick check \tab \code{nReps = 100}, serial \tab
#'     Fast; for development only \cr
#'   Standard analysis \tab \code{nReps = 1000}, parallel \tab
#'     Recommended for publication \cr
#'   Large network (\eqn{p > 15}) \tab \code{parallel = TRUE} \tab
#'     Serial runtime becomes prohibitive \cr
#'   Reproducibility \tab Always set \code{seed} \tab
#'     Report the seed value in the paper \cr
#' }
#'
#' @examples
#' ## Example 1: serial execution (small demonstration)
#' \dontrun{
#' data("single_gds")
#' gs_fit <- bootnet::estimateNetwork(single_gds,
#'             default = "IsingFit", rule = "AND")
#'
#' simResult <- stabilityNIRAtest(
#'   edge_weights               = gs_fit$graph,
#'   thresholds                 = gs_fit$intercepts,
#'   perturbation_type          = "alleviating",
#'   amount_of_SDs_perturbation = 2,
#'   nReps                      = 100,
#'   parallel                   = FALSE,
#'   seed                       = 2025
#' )
#' dim(simResult$mean)   # 100 x (p + 1)
#' }
#'
#' ## Example 2: parallel execution (recommended for publication)
#' \dontrun{
#' simResult <- stabilityNIRAtest(
#'   edge_weights               = gs_fit$graph,
#'   thresholds                 = gs_fit$intercepts,
#'   perturbation_type          = "alleviating",
#'   amount_of_SDs_perturbation = 2,
#'   nReps                      = 1000,
#'   parallel                   = TRUE,
#'   ncores                     = 6,
#'   seed                       = 2025
#' )
#' repeat_node <- findMaxN(simResult, n = 7)
#' plotMaxN(repeat_node, low = 1, high = 7)
#' }
#'
#' ## Example 3: full tutorial pipeline
#' \dontrun{
#' data("single_gds")
#' gs_fit <- bootnet::estimateNetwork(single_gds,
#'             default = "IsingFit", rule = "AND")
#'
#' # Parallel stability test
#' simResult <- stabilityNIRAtest(
#'   gs_fit$graph, gs_fit$intercepts,
#'   "alleviating", 2,
#'   nReps = 1000, parallel = TRUE, ncores = 6
#' )
#' # Compute ranking frequencies (top 7 positions)
#' repeat_node <- findMaxN(simResult, n = 7)
#' print(repeat_node)
#' # Visualise stability
#' plotMaxN(repeat_node, low = 1, high = 7)
#' }
#'
#' @author
#' Fei Wang \email{bjnwangfei0501@@outlook.com}
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{findMaxN}} — pass the returned list here to compute
#'     node ranking frequencies across repetitions.
#'   \item \code{\link{plotMaxN}} — visualise the \code{findMaxN} output
#'     as a stacked bar plot.
#'   \item \code{\link{runMgmmAnalysis}} — run this \emph{before}
#'     \code{stabilityNIRAtest} to confirm NIRA's theoretical assumptions.
#'   \item \code{\link{permutationNIRAtest}} — formal significance test
#'     to run alongside the stability test.
#'   \item \code{\link[parallel]{makeCluster}} — underlying cluster
#'     creation used in parallel mode.
#'   \item \code{\link[nodeIdentifyR]{simulateResponses}} — the simulation
#'     engine called in each repetition.
#' }
#'
#' @importFrom parallel detectCores makeCluster stopCluster
#'   clusterSetRNGStream clusterEvalQ clusterExport parLapply
#' @importFrom progress progress_bar
#' @importFrom stats sd
#' @importFrom nodeIdentifyR simulateResponses calculateSumScores
#'   prepareDFforPlottingAndANOVA
#' @export
stabilityNIRAtest <- function(edge_weights,
                               thresholds,
                               perturbation_type,
                               amount_of_SDs_perturbation,
                               nReps    = 1000,
                               parallel = FALSE,
                               ncores   = NULL,
                               seed     = 2025) {

  if (!is.null(seed)) set.seed(seed)

  simResult <- list(
    SimSamples = vector("list", nReps),
    mean       = list(),
    sd         = list()
  )

  one_rep <- function(i) {
    gs  <- nodeIdentifyR::simulateResponses(
             edge_weights, thresholds,
             perturbation_type, amount_of_SDs_perturbation)
    gss <- nodeIdentifyR::calculateSumScores(gs)
    gsl <- nodeIdentifyR::prepareDFforPlottingAndANOVA(gss)
    cats <- sort(unique(gsl[[2]]))
    list(
      sample = gsl,
      means  = sapply(cats, function(cat)
                 mean(gsl$sumscore[gsl[[2]] == cat], na.rm = TRUE)),
      sds    = sapply(cats, function(cat)
                 stats::sd(gsl$sumscore[gsl[[2]] == cat], na.rm = TRUE))
    )
  }

  if (parallel) {
    if (is.null(ncores))
      ncores <- max(1L, parallel::detectCores() - 1L)
    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterSetRNGStream(cl, seed)
    parallel::clusterEvalQ(cl, {
      library(nodeIdentifyR)
      library(dplyr)
    })
    parallel::clusterExport(
      cl,
      varlist = c("edge_weights", "thresholds",
                  "perturbation_type", "amount_of_SDs_perturbation"),
      envir   = environment())
    all_reps <- parallel::parLapply(cl, seq_len(nReps), one_rep)

  } else {
    pb <- if (interactive())
            progress::progress_bar$new(
              format = "  Simulating [:bar] :percent | ETA: :eta",
              total  = nReps, clear = FALSE, width = 65)
          else NULL
    all_reps <- vector("list", nReps)
    for (i in seq_len(nReps)) {
      all_reps[[i]] <- one_rep(i)
      if (!is.null(pb)) pb$tick()
    }
  }

  simResult$SimSamples <- lapply(all_reps, `[[`, "sample")
  simResult$mean <- do.call(rbind, lapply(all_reps, `[[`, "means"))
  simResult$sd   <- do.call(rbind, lapply(all_reps, `[[`, "sds"))

  if (length(simResult$SimSamples) > 0) {
    cats <- sort(unique(simResult$SimSamples[[1]][[2]]))
    colnames(simResult$mean) <- cats
    colnames(simResult$sd)   <- cats
  }
  simResult
}
