# compute_isotwas passes every method the same folds; this method
# must use them.

test_that("univariate_elasticnet uses given folds for every transcript", {
  X <- sim_geno(200, 30, 84)
  Y <- sim_expr(X, 3, 85)
  f <- .isotwas_folds(nrow(X), 5, 86)
  for (compress in c(TRUE, FALSE)) {
    withr::local_options(isotwas.glmnet_compress = compress)
    for (par in c(FALSE, TRUE)) {
      mod <- univariate_elasticnet(X, Y, diag(3), nfolds = 5, verbose = FALSE,
                                   par = par, n.cores = 2, seed = 86,
                                   folds = f)
      for (j in 1:3) {
        cv <- glmnet::cv.glmnet(X, Y[, j], foldid = f$foldid, alpha = .5,
                                keep = TRUE)
        expect_equal(mod$transcripts[[j]]$predicted,
                     unname(cv$fit.preval[, which.min(cv$cvm)]),
                     tolerance = 1e-10)
      }
    }
  }
})
