# multivariate_stacking on compressed rows must match the original rows.

withr::local_options(isotwas.glmnet_compress = TRUE)

test_that("multivariate_stacking is unchanged by compression", {
  X <- sim_geno(300, 50, 65)
  Y <- sim_expr(X, 3, 66)
  run <- function() multivariate_stacking(X, Y, Omega = diag(3),
                                          verbose = FALSE, seed = 67)
  # the non-negative least squares stacking step amplifies rounding a little
  expect_same_model(run(), uncompressed(run()), tol = 1e-8)
})
