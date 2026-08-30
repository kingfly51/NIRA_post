library(testthat)
library(NIRApost)

## analyticalBridgeNIRAtest --------------------------------------------------

# --------------------------------------------------------------------------
# 1) Exact-value check against an independent brute-force enumeration.
#    Small two-community network (2 source + 2 outcome, p = 4).  The test
#    re-enumerates the 2^4 spin microstates from scratch (no shared helper),
#    computes the outcome-community activation count, Boltzmann weights,
#    and the reweighted means, then compares with $exact.
# --------------------------------------------------------------------------
brute_force_bridge <- function(thr, W, out_idx, src_idx, delta_spin) {
  p <- length(thr)
  m  <- 0.5 * thr + 0.25 * (W %*% rep(1, p))
  Ws <- 0.25 * W
  S  <- as.matrix(expand.grid(rep(list(c(-1, 1)), p)))
  h  <- -0.5 * rowSums((S %*% Ws) * S) - S %*% m
  f  <- exp(-h - max(-h))
  nB <- (rowSums(S[, out_idx, drop = FALSE]) + length(out_idx)) / 2
  Z  <- sum(f)
  ma0 <- sum(nB * f) / Z
  m2  <- sum(nB^2 * f) / Z
  out <- c(original = ma0)
  for (k in seq_along(src_idx)) {
    i <- src_idx[k]
    w <- exp(delta_spin * S[, i])
    out[k + 1] <- sum(nB * f * w) / sum(f * w)
  }
  out
}

test_that("analyticalBridgeNIRAtest: exact values match brute-force enumeration", {
  thr <- c(0.4, -0.3, 0.2, -0.5)
  W   <- matrix(c(0, 0.8, 0.6, 0.2,
                  0.8, 0,  0.4, 0.5,
                  0.6, 0.4, 0,  0.7,
                  0.2, 0.5, 0.7, 0), 4, 4, byrow = TRUE)
  groups <- c("A", "A", "B", "B")

  # Use small synthetic data with the given groups (rows >= 3 to pass checks)
  set.seed(7)
  dat <- as.data.frame(matrix(rbinom(4 * 30, 1, 0.4), ncol = 4))
  colnames(dat) <- c("n1", "n2", "n3", "n4")

  res <- analyticalBridgeNIRAtest(
    data = dat, groups = groups,
    intervention_group = "A", outcome_group = "B",
    edge_weights = W, thresholds = thr,
    amount_of_SDs_perturbation = 2, nResample = 0
  )

  sdt  <- stats::sd(thr)
  ds   <- 2 * sdt
  ds_s <- 0.5 * ds
  bf   <- brute_force_bridge(thr, W, out_idx = 3:4, src_idx = 1:2,
                             delta_spin = ds_s)

  expect_equal(res$exact$mean_activation, unname(bf), tolerance = 1e-12)
})

test_that("analyticalBridgeNIRAtest: alleviating direction reverses the effect", {
  set.seed(7)
  dat <- as.data.frame(matrix(rbinom(4 * 30, 1, 0.4), ncol = 4))
  colnames(dat) <- c("n1", "n2", "n3", "n4")
  groups <- c("A", "A", "B", "B")
  # Fixed network with all-positive edge weights: direction of the effect is
  # then guaranteed (raising a source threshold cannot lower the outcome mean).
  thr <- c(0.4, -0.3, 0.2, -0.5)
  W <- matrix(c(0, 0.8, 0.6, 0.2,
                0.8, 0,  0.4, 0.5,
                0.6, 0.4, 0,  0.7,
                0.2, 0.5, 0.7, 0), 4, 4, byrow = TRUE)

  agg <- analyticalBridgeNIRAtest(dat, groups, "A", "B",
                                  perturbation_type = "aggravating",
                                  edge_weights = W, thresholds = thr,
                                  nResample = 0)
  all <- analyticalBridgeNIRAtest(dat, groups, "A", "B",
                                  perturbation_type = "alleviating",
                                  edge_weights = W, thresholds = thr,
                                  nResample = 0)

  # Original outcome-community mean is identical in both runs.
  expect_equal(agg$exact$mean_activation[1],
               all$exact$mean_activation[1], tolerance = 1e-12)
  # Aggravating (threshold up) increases the outcome activation, alleviating
  # (threshold down) decreases it.
  expect_true(all(agg$stat$effect > 0))
  expect_true(all(all$stat$effect < 0))
  # The two directions move the mean in opposite directions.
  expect_true(all(agg$exact$mean_activation[-1] >
                    agg$exact$mean_activation[1]))
  expect_true(all(all$exact$mean_activation[-1] <
                    all$exact$mean_activation[1]))
})

# --------------------------------------------------------------------------
# 2) Structural checks on the built-in single_gds with an arbitrary
#    two-community partition.
# --------------------------------------------------------------------------
test_that("analyticalBridgeNIRAtest: returns expected structure", {
  data("single_gds", package = "NIRApost")
  groups <- c(rep("A", 4), rep("B", 3))

  set.seed(2025)
  res <- analyticalBridgeNIRAtest(
    single_gds, groups, "A", "B",
    perturbation_type = "aggravating", nResample = 20
  )

  expect_true(is.list(res))
  expect_true(all(c("stat", "random", "exact", "resample_effects",
                    "plot", "parameters") %in% names(res)))

  pA <- sum(groups == "A")
  expect_equal(nrow(res$stat), pA)
  expect_true(all(c("mean_activation", "effect", "cohen_d", "se_effect",
                    "p", "p.adjust", "percenttop_1") %in% colnames(res$stat)))
  expect_equal(nrow(res$random), pA)
  expect_true(all(c("node", "effect", "p") %in% colnames(res$random)))
  expect_true(all(res$random$p >= 0 & res$random$p <= 1))
  expect_equal(nrow(res$exact), pA + 1)
  expect_equal(dim(res$resample_effects), c(20L, pA))
  expect_s3_class(res$plot, "ggplot")
  expect_equal(res$parameters$intervention_group, "A")
  expect_equal(res$parameters$outcome_group, "B")
  expect_equal(res$parameters$source_idx, 1:4)
  expect_equal(res$parameters$outcome_idx, 5:7)
})

# --------------------------------------------------------------------------
# 3) Degenerate consistency: outcome = ALL nodes should reproduce
#    analyticalNIRAtest's exact mean activation.
# --------------------------------------------------------------------------
test_that("analyticalBridgeNIRAtest: single-community degeneracy equals analyticalNIRAtest", {
  data("single_gds", package = "NIRApost")
  groups <- rep("A", ncol(single_gds))

  set.seed(2025)
  bridge <- analyticalBridgeNIRAtest(
    single_gds, groups, "A", "A",
    perturbation_type = "aggravating", nResample = 0
  )
  single <- analyticalNIRAtest(
    single_gds, "aggravating", nResample = 0
  )

  expect_equal(bridge$exact$mean_activation,
               single$exact$mean_activation, tolerance = 1e-12)
  expect_equal(bridge$stat$effect, single$stat$effect, tolerance = 1e-12)
})

# --------------------------------------------------------------------------
# 4) Argument validation
# --------------------------------------------------------------------------
test_that("analyticalBridgeNIRAtest: rejects bad groups / data", {
  data("single_gds", package = "NIRApost")
  df  <- as.data.frame(matrix(sample(0:2, 30, replace = TRUE), 10, 3))
  colnames(df) <- c("x", "y", "z")

  expect_error(analyticalBridgeNIRAtest(single_gds,
                 rep("A", ncol(single_gds) - 1), "A", "A", nResample = 0),
               "must have length p")
  expect_error(analyticalBridgeNIRAtest(df, c("A", "B", "C"), "A", "B",
                                        nResample = 0),
               "only binary")
  expect_error(analyticalBridgeNIRAtest(single_gds,
                 rep("A", ncol(single_gds)), "Z", "A", nResample = 0),
               "intervention_group")
  expect_error(analyticalBridgeNIRAtest(single_gds,
                 rep("A", ncol(single_gds)), "A", "Z", nResample = 0),
               "outcome_group")
  expect_error(analyticalBridgeNIRAtest(single_gds,
                 rep("A", ncol(single_gds)), "A", "A",
                 method = "xyzzy", nResample = 0),
               "Invalid method")
})

# --------------------------------------------------------------------------
# 5) resample_default: bootnet point estimate + bootnet resampling must be
#    centred on the point estimate (no CI offset); "fast" forces the fast
#    estimator even when parameters are supplied.
# --------------------------------------------------------------------------
test_that("analyticalBridgeNIRAtest: bootnet resampling centres CI on point estimate", {
  skip_if_not_installed("bootnet")
  skip_if_not_installed("readxl")
  set.seed(99)
  dat <- as.data.frame(matrix(rbinom(4 * 60, 1, 0.3), ncol = 4))
  colnames(dat) <- c("n1", "n2", "n3", "n4")
  groups <- c("A", "A", "B", "B")

  fit <- bootnet::estimateNetwork(dat, default = "IsingFit")
  res <- analyticalBridgeNIRAtest(dat, groups, "A", "B", "aggravating",
              edge_weights = fit$graph, thresholds = fit$intercepts,
              nResample = 50, seed = 1, parallel = FALSE)
  # Default resample_default = NULL + supplied params -> bootnet resampling.
  expect_equal(res$parameters$resample_default, "IsingFit")
  # CI should bracket the point estimate (no systematic offset).
  expect_true(all(res$stat$ci_lower <= res$stat$mean_activation))
  expect_true(all(res$stat$ci_upper >= res$stat$mean_activation))
  expect_true(all(res$stat$ci_effect_lower <= res$stat$effect))
  expect_true(all(res$stat$ci_effect_upper >= res$stat$effect))
})

test_that("analyticalBridgeNIRAtest: resample_default='fast' forces fast resampling", {
  skip_if_not_installed("bootnet")
  set.seed(99)
  dat <- as.data.frame(matrix(rbinom(4 * 60, 1, 0.3), ncol = 4))
  colnames(dat) <- c("n1", "n2", "n3", "n4")
  groups <- c("A", "A", "B", "B")

  fit <- bootnet::estimateNetwork(dat, default = "IsingFit")
  res <- analyticalBridgeNIRAtest(dat, groups, "A", "B", "aggravating",
              edge_weights = fit$graph, thresholds = fit$intercepts,
              nResample = 50, seed = 1, parallel = FALSE,
              resample_default = "fast")
  expect_equal(res$parameters$resample_default, "fast")
})

# --------------------------------------------------------------------------
# 6) Reproducibility with a fixed seed
# --------------------------------------------------------------------------
test_that("analyticalBridgeNIRAtest: seed gives reproducible results", {
  data("single_gds", package = "NIRApost")
  groups <- c(rep("A", 4), rep("B", 3))

  a <- analyticalBridgeNIRAtest(single_gds, groups, "A", "B",
                                "aggravating", nResample = 20, seed = 42)
  b <- analyticalBridgeNIRAtest(single_gds, groups, "A", "B",
                                "aggravating", nResample = 20, seed = 42)
  expect_identical(a$stat$effect, b$stat$effect)
  expect_identical(a$resample_effects, b$resample_effects)
})
