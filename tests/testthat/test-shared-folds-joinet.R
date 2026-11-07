# compute_isotwas passes every method the same folds; this method
# must use them.

test_that("multivariate_joinet uses given folds as its external folds", {
  skip_if_not_installed("joinet", "1.0.0")
  X <- sim_geno(200, 30, 87)
  Y <- sim_expr(X, 3, 88)
  f <- .isotwas_folds(nrow(X), 5, 89)
  utils::capture.output(ref <- withr::with_seed(89, joinet::cv.joinet(
    Y = as.data.frame(Y), X = X, foldid.ext = f$foldid, nfolds.int = 5,
    alpha.base = .5, alpha.meta = .5, cvpred = TRUE)))
  mod <- multivariate_joinet(X, Y, diag(3), nfolds = 5, verbose = FALSE,
                             seed = 89, nfolds.int = 5, folds = f)
  for (j in 1:3) {
    expect_equal(mod$transcripts[[j]]$predicted, ref$cvpred[, j],
                 tolerance = 1e-10)
  }
})
