test_that("permutation input is independent of column order and extra columns", {
  x <- data.frame(sample = rep(c("original", "A", "B"), each = 4),
                  sumscore = c(0:3, 2:5, 4:7))
  a <- permutationNIRAtest(x)
  x$id <- seq_len(nrow(x))
  b <- permutationNIRAtest(x[, c("sumscore", "id", "sample")])
  expect_equal(a, b)
  expect_equal(a$stat$effect, c(2, 4))
  expect_error(permutationNIRAtest(x[x$sample != "original", ]), "original")
})

test_that("ranking handles one requested position and one repetition", {
  m <- matrix(c(0, 0, 1, 2, 2, 1), 2,
              dimnames = list(NULL, c("original", "A", "B")))
  expect_equal(findMaxN(list(mean = m), n = 1)$percenttop_1, c(.5, .5))
  expect_equal(findMaxN(list(mean = m[1, , drop = FALSE]), n = 1)$percenttop_1, c(0, 1))
})

test_that("imputation preserves numeric data and forwards explicit tuning arguments", {
  x <- data.frame(A = c(0, 1, 1, NA), B = 1:4)
  expect_type(imputeData(x, method = "mode")$A, "double")
  expect_false(anyNA(imputeData(as.matrix(x), method = "mode")))
  expect_false(anyNA(imputeData(x, method = "knn", k = 2)))
  set.seed(42)
  y <- data.frame(A = rnorm(60), B = rnorm(60), C = rnorm(60))
  y$A[c(1, 5, 10)] <- NA
  expect_false(anyNA(imputeData(y, method = "mi", m = 2, maxit = 1)))
  expect_false(anyNA(imputeData(y, method = "random_forest", ntree = 5, maxiter = 1)))
})

test_that("analytical effects returned to users agree with reported uncertainty", {
  set.seed(93)
  x <- matrix(rbinom(1200, 1, .4), ncol = 3)
  a <- analyticalNIRAtest(x, nResample = 20)
  b <- analyticalBridgeNIRAtest(x, groups = c("A", "A", "B"),
                              intervention_group = "A", outcome_group = "B", nResample = 20)
  for (r in list(a, b)) {
    expect_equal(unname(apply(r$resample_effects, 2, sd)), r$stat$se_effect)
    expect_equal(unname(apply(r$resample_effects, 2, quantile, .025)), r$stat$ci_effect_lower)
    expect_equal(sum(r$stat$percenttop_1), 1)
    expect_false("random" %in% names(r))
  }
})

test_that("large finite perturbations have finite independent-node expectations", {
  x <- as.matrix(expand.grid(A = 0:1, B = 0:1))
  r <- analyticalNIRAtest(x, thresholds = c(-1, 1), edge_weights = matrix(0, 2, 2),
                         amount_of_SDs_perturbation = 2000, nResample = 0)
  expect_equal(r$exact$mean_activation, c(1, 1 + plogis(1), 1 + plogis(-1)))
  r <- analyticalBridgeNIRAtest(x, c("A", "B"), "A", "B", thresholds = c(-1, 1),
                                edge_weights = matrix(0, 2, 2),
                                amount_of_SDs_perturbation = 2000, nResample = 0)
  expect_equal(r$stat$effect, 0, tolerance = 1e-12)
  expect_error(analyticalNIRAtest(cbind(A = rep(0, 8), B = rep(0:1, 4)), nResample = 0), "finite")
})
