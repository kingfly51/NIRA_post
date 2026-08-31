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

## permutationNIRAtest --------------------------------------------------------
test_that("permutationNIRAtest: returns a 'random' comparison element", {
  set.seed(1)
  # 列顺序与 prepareDFforPlottingAndANOVA 一致：[sumscore, sample]
  df <- data.frame(
    sumscore = c(rnorm(40, 10, 2), rnorm(40, 12, 2),
                 rnorm(40, 11, 2), rnorm(40, 10.5, 2)),
    sample   = c(rep("original", 40), rep("A", 40),
                 rep("B", 40), rep("C", 40))
  )
  r <- permutationNIRAtest(df, method = "bonferroni")
  expect_true("random" %in% names(r))
  expect_true(all(c("node", "effect", "p") %in% colnames(r$random)))
  expect_true(all(r$random$p >= 0 & r$random$p <= 1))
  # 效应最大的节点 p 值最小（最好靶点）
  expect_equal(r$random$node[which.max(r$random$effect)],
               r$random$node[which.min(r$random$p)])
})
test_that("permutationNIRAtest: direction aligns random test with intervention direction", {
  set.seed(1)
  # A 是最有效的 alleviating 靶点（最负效应），C 为反向（正效应）
  df <- data.frame(
    sumscore = c(rnorm(60, 10, 2), rnorm(60, 6, 2),
                 rnorm(60, 8, 2), rnorm(60, 11, 2)),
    sample   = c(rep("original", 60), rep("A", 60),
                 rep("B", 60), rep("C", 60))
  )
  r_abs <- permutationNIRAtest(df, direction = "abs")
  r_all <- permutationNIRAtest(df, direction = "alleviating")
  # 默认 abs 向后兼容
  expect_equal(r_abs$random$p, r_all$random$p) # C 的反向效应在两种记分下都非最优
  # 方向对齐下，最负效应节点（A）p 最小
  expect_equal(r_all$random$node[which.min(r_all$random$p)], "A")
  # 非法 direction 报错
  expect_error(permutationNIRAtest(df, direction = "bogus"))
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

## runMgmmAnalysis ------------------------------------------------------------
test_that("runMgmmAnalysis: accepts ncores parameter", {
  expect_true("ncores" %in% names(formals(runMgmmAnalysis)))
})
test_that("runMgmmAnalysis: rejects invalid ncores", {
  expect_error(runMgmmAnalysis(matrix(0, 4, 3), ncores = 0))
  expect_error(runMgmmAnalysis(matrix(0, 4, 3), ncores = 1.5))
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
