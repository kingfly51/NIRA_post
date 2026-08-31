library(testthat)
library(NIRApost)

## analyticalNIRAtest ---------------------------------------------------------
data("single_gds", package = "NIRApost")

test_that("analyticalNIRAtest: exact mean activation matches Isinglandr", {
  skip_if_not_installed("Isinglandr")
  set.seed(2025)
  r <- analyticalNIRAtest(single_gds, "aggravating", nResample = 0)
  thr <- r$parameters$thresholds
  W   <- r$parameters$edge_weights

  ma_isl <- function(tt, Ww) {
    l <- Isinglandr::make_2d_Isingland(thresholds = tt, weiadj = Ww,
                                       transform = TRUE)
    sum(l$dist$n_active * l$dist$sum_freq) / sum(l$dist$sum_freq)
  }
  for (i in seq_len(nrow(r$exact))) {
    cond <- r$exact$condition[i]
    if (cond == "original") {
      f <- ma_isl(thr, W)
    } else {
      tt <- thr
      tt[which(colnames(single_gds) == cond)] <-
        tt[which(colnames(single_gds) == cond)] + 2 * stats::sd(thr)
      f <- ma_isl(tt, W)
    }
    expect_equal(r$exact$mean_activation[i], f, tolerance = 1e-10)
  }
})

test_that("analyticalNIRAtest: returns expected structure", {
  set.seed(2025)
  r <- analyticalNIRAtest(single_gds, "aggravating", nResample = 20)
  expect_true(is.list(r))
  expect_true(all(c("stat", "random", "exact", "resample_effects",
                    "plot", "parameters") %in% names(r)))
  p <- ncol(single_gds)
  expect_equal(nrow(r$stat), p)
  expect_true(all(c("mean_activation", "effect", "se_effect", "p",
                    "p.adjust", "percenttop_1") %in% colnames(r$stat)))
  expect_equal(nrow(r$random), p)
  expect_true(all(c("node", "effect", "p") %in% colnames(r$random)))
  expect_true(all(r$random$p >= 0 & r$random$p <= 1))
  expect_equal(nrow(r$exact), p + 1)
  expect_equal(dim(r$resample_effects), c(20L, p))
  expect_s3_class(r$plot, "ggplot")
})

test_that("analyticalNIRAtest: random-target p-value ranks the best node lowest", {
  set.seed(2025)
  r <- analyticalNIRAtest(single_gds, "aggravating", nResample = 20)
  best <- r$random$node[which.max(r$random$effect)]
  expect_equal(best, "ASH")
  expect_equal(r$random$p[r$random$node == best],
               min(r$random$p))
})

test_that("analyticalNIRAtest: random-target p-values have discriminating power", {
  # The random-target test must NOT collapse to p = 1 for every node (the
  # symptom of comparing a point estimate against a resampled pool).  Here we
  # verify the p-value range spans well below 1 and ranks by effect size.
  set.seed(2025)
  r <- analyticalNIRAtest(single_gds, "aggravating", nResample = 50)
  pvals <- r$random$p
  expect_true(min(pvals) < 0.1, info = paste("min p =", min(pvals)))
  expect_true(max(pvals) > 0.5, info = paste("max p =", max(pvals)))
  # Rank order: larger (signed) effect -> smaller p
  ord <- order(r$random$effect, decreasing = TRUE)
  expect_true(all(diff(r$random$p[ord]) >= 0) ||
                # tolerate tiny ties from discrete counts
                all(diff(r$random$p[ord]) >= -0.001))
})

test_that("analyticalNIRAtest: ASH is the strongest aggravating target", {
  set.seed(2025)
  r <- analyticalNIRAtest(single_gds, "aggravating", nResample = 20)
  expect_equal(r$stat["ASH", "effect"], max(r$stat$effect))
})

test_that("analyticalNIRAtest: nResample = 0 gives exact-only output", {
  r <- analyticalNIRAtest(single_gds, "aggravating", nResample = 0)
  expect_null(r$resample_effects)
  expect_true(all(is.na(r$stat$p)))
  expect_true(all(is.na(r$stat$percenttop_1)))
  expect_equal(nrow(r$exact), ncol(single_gds) + 1)
})

test_that("analyticalNIRAtest: seed gives reproducible results", {
  a <- analyticalNIRAtest(single_gds, "aggravating", nResample = 20, seed = 42)
  b <- analyticalNIRAtest(single_gds, "aggravating", nResample = 20, seed = 42)
  expect_identical(a$stat$effect, b$stat$effect)
  expect_identical(a$resample_effects, b$resample_effects)
})

test_that("analyticalNIRAtest: rejects non-binary data", {
  df <- as.data.frame(matrix(sample(0:2, 30, replace = TRUE), 10, 3))
  expect_error(analyticalNIRAtest(df, "aggravating", nResample = 0))
})

test_that("analyticalNIRAtest: requires both or neither of parameters", {
  W <- matrix(0, 3, 3)
  expect_error(analyticalNIRAtest(single_gds[, 1:3], "aggravating",
                                  edge_weights = W, nResample = 0),
               "BOTH 'thresholds' and 'edge_weights'")
  expect_error(analyticalNIRAtest(single_gds[, 1:3], "aggravating",
                                  thresholds = c(-1, -1, -1), nResample = 0),
               "BOTH 'thresholds' and 'edge_weights'")
})

test_that("analyticalNIRAtest: rejects invalid method", {
  expect_error(analyticalNIRAtest(single_gds, "aggravating",
                                  method = "xyzzy", nResample = 0))
})

test_that("analyticalNIRAtest: bootnet resampling centres CI on point estimate", {
  skip_if_not_installed("bootnet")
  set.seed(99)
  dat <- as.data.frame(matrix(rbinom(4 * 60, 1, 0.3), ncol = 4))
  colnames(dat) <- c("n1", "n2", "n3", "n4")

  fit <- bootnet::estimateNetwork(dat, default = "IsingFit")
  res <- analyticalNIRAtest(dat, "aggravating",
              edge_weights = fit$graph, thresholds = fit$intercepts,
              nResample = 50, seed = 1, parallel = FALSE)
  # Default resample_default = NULL + supplied params -> bootnet resampling.
  expect_equal(res$parameters$resample_default, "IsingFit")
  # CI should bracket the point estimate (no systematic offset).
  expect_true(all(res$stat$ci_lower <= res$stat$mean_activation, na.rm = TRUE))
  expect_true(all(res$stat$ci_upper >= res$stat$mean_activation, na.rm = TRUE))
  expect_true(all(res$stat$ci_effect_lower <= res$stat$effect, na.rm = TRUE))
  expect_true(all(res$stat$ci_effect_upper >= res$stat$effect, na.rm = TRUE))
})

test_that("analyticalNIRAtest: resample_default='fast' forces fast resampling", {
  skip_if_not_installed("bootnet")
  set.seed(99)
  dat <- as.data.frame(matrix(rbinom(4 * 60, 1, 0.3), ncol = 4))
  colnames(dat) <- c("n1", "n2", "n3", "n4")

  fit <- bootnet::estimateNetwork(dat, default = "IsingFit")
  res <- analyticalNIRAtest(dat, "aggravating",
              edge_weights = fit$graph, thresholds = fit$intercepts,
              nResample = 50, seed = 1, parallel = FALSE,
              resample_default = "fast")
  expect_equal(res$parameters$resample_default, "fast")
})
