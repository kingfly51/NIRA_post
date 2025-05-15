#' @param original To simulate the original sample group, a column of vectors should be input
#' @param group For the sample group after simulation intervention, a column of vectors should be input
#' @param n_perm Number of permutations,5000
#' @examples permutation_test(original, group, 5000)
#' @author Written by Wang Fei 2025/03/22
#' NIRA
#' State Key Laboratory of Cognitive Neuroscience and Learning,Beijing Normal University
#' @email bjnwangfei0501@outlook.com

permutation_test <- function(original, group, n_perm) {
  
  stopifnot(n_perm > 0)
  
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
