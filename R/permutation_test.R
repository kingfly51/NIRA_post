#' Permutation Test for Group Differences
#'
#' @description
#' Performs a permutation test to compare means between two groups.
#'
#' @param original Numeric vector of values from original sample.
#' @param group Numeric vector of values from intervention sample.
#' @param n_perm Number of permutations (default: 5000).
#'
#' @return
#' A list containing:
#' - `mean_original`, `mean_other`: Group means
#' - `sd_*`, `se_*`: Standard deviations and standard errors
#' - `ci_*`, `ciLower_*`, `ciUpper_*`: 95% confidence intervals
#' - `pooled_std_dev`: Pooled standard deviation
#' - `cohens_d`: Cohen's d effect size
#' - `p_value`: Permutation p-value
#' - `n_permutations`: Number of permutations
#' - `permutation_distribution`: All permuted differences
#'
#' @examples
#' orig <- rnorm(100, mean = 5)
#' interv <- rnorm(100, mean = 5.5)
#' results <- permutation_test(orig, interv, 1000)
#'
#' @author
#' Wang Fei, State Key Laboratory of Cognitive Neuroscience and Learning,
#' Beijing Normal University (2025-03-22)
#'
#' @importFrom stats var sd
#' @export

permutation_test <- function(original, group, n_perm = 5000) {
  # --- Input validation ---
  if (!is.numeric(original) || !is.numeric(group)) {
    stop("Both inputs must be numeric vectors")
  }

  if (length(original) < 3 || length(group) < 3) {
    stop("Each group should have at least 3 observations")
  }

  if (!is.numeric(n_perm) || n_perm < 100) {
    stop("Number of permutations should be ≥ 100")
  }

  n1 <- length(original)
  n2 <- length(group)

  mean_original <- mean(original)
  mean_group <- mean(group)
  observed_diff <- mean_original - mean_group

  sd_original <- sd(original)
  sd_group <- sd(group)

  se_original <- sd_original/sqrt(n1)
  se_group <- sd_group/sqrt(n2)

  ci_original <- 1.96*se_original
  ci_group <- 1.96*se_group

  ciLower_original <- mean_original - ci_original
  ciUpper_original <- mean_original + ci_original

  ciLower_group <- mean_group - ci_group
  ciUpper_group <- mean_group + ci_group

  var_original <- var(original)
  var_group <- var(group)

  s_pooled <- sqrt(((n1 - 1) * var_original + (n2 - 1) * var_group) / (n1 + n2 - 2))
  cohens_d <- observed_diff / s_pooled

  combined <- c(original, group)
  total_n <- n1 + n2
  perm_diffs <- numeric(n_perm)

  for (i in 1:n_perm) {
    shuffled <- sample(combined, total_n, replace = FALSE)
    new_original <- shuffled[1:n1]
    new_group <- shuffled[(n1 + 1):total_n]
    perm_diffs[i] <- mean(new_original) - mean(new_group)
  }

  extreme_count <- sum(abs(perm_diffs) >= abs(observed_diff))
  p_value <- (extreme_count + 1) / (n_perm + 1)

  permu_list <- list(
    mean_original = mean_original,
    mean_other = mean_group,
    sd_original = sd_original,
    sd_other = sd_group,
    se_original = se_original,
    se_other = se_group,
    ci_original = ci_original,
    ci_other = ci_group,
    ciLower_original = ciLower_original,
    ciUpper_original = ciUpper_original,
    ciLower_other = ciLower_group,
    ciUpper_other = ciUpper_group,
    pooled_std_dev = s_pooled,
    cohens_d = cohens_d,
    p_value = p_value,
    n_permutations = n_perm,
    permutation_distribution = perm_diffs
  )
  return(permu_list)
}
