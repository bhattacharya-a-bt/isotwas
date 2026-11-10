# compute_isotwas passes every method the same folds; this method
# must use them.

test_that("multivariate_mtlasso scores R2 on the given folds", {
  X <- sim_geno(150, 30, 93)
  Y <- sim_expr(X, 3, 94)
  f <- .isotwas_folds(nrow(X), 5, 95)
  run <- function(folds) multivariate_mtlasso(X, Y, nlambda = 5, nfolds = 5,
                                              seed = 95, folds = folds)
  a <- run(f)
  # the R2 folds used to be a second draw after the lambda folds
  d <- run(NULL)
  expect_identical(a$transcripts[[1]]$weights, d$transcripts[[1]]$weights)
  expect_false(identical(a$transcripts[[1]]$predicted,
                         d$transcripts[[1]]$predicted))
  expect_identical(a$best_lambda, d$best_lambda)
  # with a fixed lambda only the R2 folds are drawn; given folds replace
  # them, in either form
  b1 <- multivariate_mtlasso(X, Y, lambda = a$best_lambda, seed = 95,
                             folds = f)
  b2 <- multivariate_mtlasso(X, Y, lambda = a$best_lambda, seed = 95,
                             folds = f$foldid)
  expect_identical(b1$transcripts[[1]]$predicted, a$transcripts[[1]]$predicted)
  expect_identical(b2$transcripts[[1]]$predicted, a$transcripts[[1]]$predicted)
})
