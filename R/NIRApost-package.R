#' NIRApost: Post-Processing Tools for the NodeIdentifyR Algorithm (NIRA)
#'
#' @description
#' Three-step validation pipeline for NIRA in Ising network analysis:
#' \enumerate{
#'   \item Moderation effect testing (prerequisite):
#'     \code{\link{runMgmmAnalysis}}
#'   \item Permutation significance testing:
#'     \code{\link{permutationTest}}, \code{\link{permutationNIRAtest}},
#'     \code{\link{plotMeanNIRA}}
#'   \item Repeated-simulation stability assessment:
#'     \code{\link{stabilityNIRAtest}}, \code{\link{findMaxN}},
#'     \code{\link{plotMaxN}}
#' }
#'
#' @section Data preprocessing utilities:
#' \code{\link{checkMissing}}, \code{\link{imputeData}}
#'
#' @section Datasets:
#' \code{single_gds} (GAD-7, N=2404), \code{PHQ9} (N=45829),
#' \code{GDS9} (N=3097), \code{ACE} (N=525), \code{disease} (N=27353)
#'
#' @author Fei Wang \email{bjnwangfei0501@@outlook.com},
#'   Yiming Wu, Yibo Wu, Tingshao Zhu
#' @keywords internal
"_PACKAGE"

## Suppress R CMD check NOTE: no visible binding for global variable
utils::globalVariables(c(
  "sample", "sumscore", "row",
  "node", "nt", "mean", "newid", "ciLower", "ciUpper", "p.adjust", "stars",
  "Pct", "Rank"
))

