library(testthat)
library(NIRApost)

## checkMissing ---------------------------------------------------------------
test_that("checkMissing: identifies complete and incomplete columns", {
  df <- data.frame(A = 1:5,
                   B = c(NA, 2, NA, 4, 5),
                   C = c(1, 2, 3, NA, NA))
  r  <- checkMissing(df)
  expect_equal(r$complete_cols, "A")
  expect_setequal(r$incomplete_cols, c("B", "C"))
  expect_true("mi"  %in% r$available_methods)
  expect_true("knn" %in% r$available_methods)
})
test_that("checkMissing: restricts methods when all cols have NAs", {
  r <- checkMissing(data.frame(A = c(1, NA, 3), B = c(NA, 2, 3)))
  expect_length(r$complete_cols, 0)
  expect_setequal(r$available_methods, c("mode", "random_forest", "mi"))
})
test_that("checkMissing: rejects non-data-frame input", {
  expect_error(checkMissing("not a data frame"))
})

## permutationTest ------------------------------------------------------------
test_that("permutationTest: returns correct list structure", {
  set.seed(42)
  r <- permutationTest(rnorm(60, 5), rnorm(60, 6), nPerm = 300)
  expect_true(is.list(r))
  expect_equal(r$n_permutations, 300)
  expect_gte(r$p_value, 0)
  expect_lte(r$p_value, 1)
  expect_true("cohens_d" %in% names(r))
  expect_length(r$permutation_distribution, 300)
})
test_that("permutationTest: rejects non-numeric input", {
  expect_error(permutationTest(letters[1:10], 1:10))
})
test_that("permutationTest: rejects groups with fewer than 3 obs", {
  expect_error(permutationTest(c(1, 2), c(1, 2, 3)))
})
test_that("permutationTest: rejects nPerm < 100", {
  expect_error(permutationTest(rnorm(20), rnorm(20), nPerm = 50))
})

## findMaxN -------------------------------------------------------------------
test_that("findMaxN: returns correct dimensions and valid proportions", {
  set.seed(1)
  node_names <- paste0("node", 1:5)
  mock_sim <- list(
    SimSamples = vector("list", 10),
    mean = matrix(rnorm(10 * 6), nrow = 10,
                  dimnames = list(NULL, c("original", node_names))),
    sd   = matrix(abs(rnorm(10 * 6)), nrow = 10,
                  dimnames = list(NULL, c("original", node_names)))
  )
  r <- findMaxN(mock_sim, n = 3)
  expect_equal(nrow(r), 5)
  expect_true("percenttop_1" %in% colnames(r))
  expect_true("repeattop_3"  %in% colnames(r))
  pct_cols <- grep("percent", colnames(r), value = TRUE)
  expect_true(all(r[, pct_cols] >= 0))
  expect_true(all(r[, pct_cols] <= 1))
})
test_that("findMaxN: rejects missing original column", {
  bad_sim <- list(mean = matrix(1:12, 3, 4,
                   dimnames = list(NULL, paste0("n", 1:4))))
  expect_error(findMaxN(bad_sim, n = 2))
})

## stabilityNIRAtest ----------------------------------------------------------
test_that("stabilityNIRAtest: accepts nReps parameter", {
  expect_true("nReps" %in% names(formals(stabilityNIRAtest)))
})

## imputeData -----------------------------------------------------------------
test_that("imputeData: returns original data when no NAs", {
  df <- data.frame(A = 1:5, B = 6:10)
  expect_equal(imputeData(df, method = "mode"), df)
})
test_that("imputeData: rejects invalid method", {
  expect_error(imputeData(data.frame(A = c(1, NA, 1)), method = "xyzzy"))
})
test_that("imputeData: rejects non-data-frame input", {
  expect_error(imputeData("not a data frame"))
})
