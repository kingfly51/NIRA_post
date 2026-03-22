#' Two-Sample Permutation Test for Differences in Means
#'
#' @description
#' Tests whether the means of two groups differ significantly by randomly
#' reassigning observations between groups many times and comparing the
#' observed mean difference against the resulting null distribution.
#' Unlike parametric \emph{t}-tests, no distributional assumptions are
#' required, making this approach well-suited to the simulated sum-score
#' data produced by the NodeIdentifyR Algorithm (NIRA).
#'
#' This function is called internally by \code{\link{permutationNIRAtest}}
#' for every node in the network, but it can also be used as a
#' stand-alone test for any two numeric vectors.
#'
#' The p-value uses the \emph{plus-one} convention
#' (\eqn{p = (k + 1) / (B + 1)}, where \eqn{k} is the number of
#' permuted differences at least as extreme as the observed difference
#' and \eqn{B} is the number of permutations) to avoid a zero p-value
#' and control the Type-I error rate.
#'
#' @param original Numeric vector of scores from the \strong{original}
#'   (unperturbed) simulated network.  Must contain at least 3 observations.
#' @param group Numeric vector of scores from the \strong{intervention}
#'   simulated network (one node perturbed).  Must contain at least 3
#'   observations.
#' @param nPerm Positive integer specifying the number of permutations.
#'   Must be \eqn{\geq 100}.  Default: \code{5000}.  Use at least
#'   \code{5000} for publication-quality p-values; \code{1000} is
#'   sufficient for exploratory screening.
#'
#' @return
#' A named \code{list} with 17 elements:
#'
#' \describe{
#'   \item{\code{mean_original}}{Mean of the original group.}
#'   \item{\code{mean_other}}{Mean of the intervention group.}
#'   \item{\code{sd_original}}{Standard deviation of the original group.}
#'   \item{\code{sd_other}}{Standard deviation of the intervention group.}
#'   \item{\code{se_original}}{Standard error of the original group mean.}
#'   \item{\code{se_other}}{Standard error of the intervention group mean.}
#'   \item{\code{ci_original}}{Half-width of the 95 \% CI for the original
#'     mean (\eqn{1.96 \times \mathrm{SE}}).}
#'   \item{\code{ci_other}}{Half-width of the 95 \% CI for the
#'     intervention mean.}
#'   \item{\code{ciLower_original}}{Lower 95 \% CI bound for the original
#'     mean.}
#'   \item{\code{ciUpper_original}}{Upper 95 \% CI bound for the original
#'     mean.}
#'   \item{\code{ciLower_other}}{Lower 95 \% CI bound for the intervention
#'     mean.}
#'   \item{\code{ciUpper_other}}{Upper 95 \% CI bound for the intervention
#'     mean.}
#'   \item{\code{pooled_std_dev}}{Pooled standard deviation used to
#'     compute Cohen's \eqn{d}.}
#'   \item{\code{cohens_d}}{Cohen's \eqn{d} effect size
#'     (\eqn{(\bar{x}_\text{orig} - \bar{x}_\text{other}) / s_p}).}
#'   \item{\code{p_value}}{Two-sided permutation p-value.}
#'   \item{\code{n_permutations}}{Number of permutations actually used.}
#'   \item{\code{permutation_distribution}}{Numeric vector of length
#'     \code{nPerm} containing all permuted mean differences, useful for
#'     diagnostic plotting.}
#' }
#'
#' @section Method Recommendations:
#' \tabular{lll}{
#'   \strong{Purpose} \tab \strong{nPerm} \tab \strong{Note} \cr
#'   Exploratory / screening \tab 1 000 \tab Fast; acceptable precision \cr
#'   Standard analysis \tab 5 000 \tab Default; good balance \cr
#'   Publication \tab 10 000+ \tab Most stable p-values \cr
#' }
#'
#' @examples
#' ## Example 1: detect a real group difference
#' set.seed(42)
#' result <- permutationTest(
#'   original = rnorm(200, mean = 5, sd = 1),
#'   group    = rnorm(200, mean = 6, sd = 1),
#'   nPerm    = 1000
#' )
#' result$p_value     # should be < 0.05
#' result$cohens_d    # effect size
#'
#' ## Example 2: no group difference (null scenario)
#' set.seed(1)
#' result_null <- permutationTest(
#'   original = rnorm(100, mean = 5),
#'   group    = rnorm(100, mean = 5),
#'   nPerm    = 1000
#' )
#' result_null$p_value   # should be large
#'
#' ## Example 3: inspect the permutation null distribution
#' \dontrun{
#' set.seed(7)
#' res <- permutationTest(rnorm(150, 4), rnorm(150, 5), nPerm = 2000)
#' hist(res$permutation_distribution,
#'      main = "Permutation null distribution",
#'      xlab = "Mean difference")
#' abline(v = res$mean_original - res$mean_other, col = "red", lwd = 2)
#' }
#'
#' @author
#' Fei Wang \email{bjnwangfei0501@@outlook.com}
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{permutationNIRAtest}} — applies this function to all
#'     nodes in a NIRA result and corrects for multiple comparisons.
#'   \item \code{\link{plotMeanNIRA}} — visualises the output of
#'     \code{permutationNIRAtest} with significance stars derived from this
#'     function's p-values.
#'   \item \code{\link[stats]{p.adjust}} — multiple-comparison correction
#'     applied to the p-values returned here.
#' }
#'
#' @importFrom stats var sd
#' @export
permutationTest <- function(original, group, nPerm = 5000) {

  if (!is.numeric(original) || !is.numeric(group))
    stop("Both 'original' and 'group' must be numeric vectors.")
  if (length(original) < 3 || length(group) < 3)
    stop("Each group must contain at least 3 observations.")
  if (!is.numeric(nPerm) || nPerm < 100)
    stop("nPerm must be a numeric value >= 100.")

  n1 <- length(original)
  n2 <- length(group)

  mo  <- mean(original);  mg  <- mean(group)
  obs <- mo - mg

  so  <- stats::sd(original);  sg  <- stats::sd(group)
  seo <- so / sqrt(n1);        seg <- sg / sqrt(n2)

  vo  <- stats::var(original); vg  <- stats::var(group)
  sp  <- as.numeric(sqrt(((n1 - 1) * vo + (n2 - 1) * vg) / (n1 + n2 - 2)))

  combined  <- c(original, group)
  perm_diff <- numeric(nPerm)
  for (i in seq_len(nPerm)) {
    lb          <- sample(rep(c(1L, 2L), times = c(n1, n2)))
    perm_diff[i] <- mean(combined[lb == 1L]) - mean(combined[lb == 2L])
  }
  extreme_count <- sum(abs(perm_diff) >= abs(obs))

  list(
    mean_original            = mo,
    mean_other               = mg,
    sd_original              = so,
    sd_other                 = sg,
    se_original              = seo,
    se_other                 = seg,
    ci_original              = 1.96 * seo,
    ci_other                 = 1.96 * seg,
    ciLower_original         = mo - 1.96 * seo,
    ciUpper_original         = mo + 1.96 * seo,
    ciLower_other            = mg - 1.96 * seg,
    ciUpper_other            = mg + 1.96 * seg,
    pooled_std_dev           = sp,
    cohens_d                 = obs / sp,
    p_value                  = (extreme_count + 1) / (nPerm + 1),
    n_permutations           = nPerm,
    permutation_distribution = perm_diff
  )
}