# .cv_joinet() and .joinet_fit() must reproduce joinet 1.0.0's cv.joinet()
# and joinet() for gaussian outputs, random draws included.

withr::local_options(isotwas.glmnet_compress = TRUE)

test_that(".cv_joinet and .joinet_fit match joinet 1.0.0", {
  skip_if_not_installed("joinet", "1.0.0")
  X <- sim_geno(250, 40, 61)
  Y <- sim_expr(X, 3, 62)
  ref_pred <- withr::with_seed(63, suppressWarnings(utils::capture.output(
    r <- joinet::cv.joinet(Y = as.data.frame(Y), X = X, nfolds.ext = 5,
                           nfolds.int = 5, alpha.base = .5, alpha.meta = .5,
                           cvpred = TRUE))))
  new_pred <- withr::with_seed(63, .cv_joinet(X, Y, 5, 5, .5, .5))
  expect_equal(new_pred, r$cvpred, tolerance = 1e-10)

  ref <- withr::with_seed(64, joinet::joinet(Y = Y, X = X, nfolds = 5,
                                             alpha.base = .5, alpha.meta = .5))
  new <- withr::with_seed(64, .joinet_fit(X, Y, .joinet_folds(Y[, 1], 5),
                                          .5, .5))
  expect_equal(unname(.joinet_coef(new)$beta),
               unname(coef(ref)$beta), tolerance = 1e-10)
  expect_equal(unname(.joinet_coef(new)$alpha),
               unname(coef(ref)$alpha), tolerance = 1e-10)
  expect_equal(.joinet_predict(new, X[1:10, ]),
               predict(ref, X[1:10, ]), tolerance = 1e-10)
})

test_that("multivariate_joinet is unchanged by compression", {
  X <- sim_geno(250, 40, 57)
  Y <- sim_expr(X, 3, 58)
  run <- function() multivariate_joinet(X, Y, diag(3), nfolds = 5,
                                        verbose = FALSE, seed = 59,
                                        nfolds.int = 5)
  expect_same_model(run(), uncompressed(run()))
})
