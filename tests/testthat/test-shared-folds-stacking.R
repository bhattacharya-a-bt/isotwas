# compute_isotwas passes every method the same folds; this method
# must use them.

test_that("multivariate_stacking uses given folds for the stacking weights", {
  X <- sim_geno(200, 30, 90)
  Y <- sim_expr(X, 3, 91)
  f <- .isotwas_folds(nrow(X), 5, 92)
  run <- function(folds) multivariate_stacking(X, Y, Omega = diag(3),
                                               seed = 92, folds = folds)
  a <- run(f)
  b <- run(f$foldid)
  expect_identical(a$stack_weights_matrix, b$stack_weights_matrix)
  expect_identical(a$transcripts[[1]]$predicted, b$transcripts[[1]]$predicted)
  # the folds stacking drew itself came after the tuning draws, so differ
  expect_false(identical(a$transcripts[[1]]$predicted,
                         run(NULL)$transcripts[[1]]$predicted))
})
