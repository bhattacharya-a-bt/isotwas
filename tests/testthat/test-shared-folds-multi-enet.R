# compute_isotwas passes every method the same folds; this method
# must use them.

test_that("multivariate_elasticnet uses given folds for all three alphas", {
  X <- sim_geno(200, 30, 81)
  Y <- sim_expr(X, 3, 82)
  f <- .isotwas_folds(nrow(X), 5, 83)
  mod <- multivariate_elasticnet(X, Y, diag(3), nfolds = 5, verbose = FALSE,
                                 seed = 83, folds = f)
  ref <- lapply(c(0, .5, 1), function(a) {
    cv <- glmnet::cv.glmnet(X, Y, foldid = f$foldid, family = "mgaussian",
                            alpha = a, standardize = FALSE, keep = TRUE)
    cv$fit.preval[, , which.min(cv$cvm)]
  })
  for (j in 1:3) {
    d <- sapply(ref, function(r) max(abs(r[, j] - mod$transcripts[[j]]$predicted)))
    expect_lt(min(d), 1e-10)
  }
  lab <- multivariate_elasticnet(X, Y, diag(3), nfolds = 5, verbose = FALSE,
                                 seed = 83, folds = f$foldid)
  expect_identical(lab$transcripts[[1]]$predicted,
                   mod$transcripts[[1]]$predicted)
})
