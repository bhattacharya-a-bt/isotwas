# multivariate_sgl's compiled solver must reproduce the R solver it replaces.

test_that("the compiled sparse group lasso matches .fit_sgl_r", {
  X <- scale(sim_geno(120, 60, 1), scale = FALSE)
  Y <- scale(sim_expr(X, 3, 2), scale = FALSE)
  cd_gram <- .cd_data(X, Y, NULL, crossprod(X), crossprod(X, Y))
  B0 <- withr::with_seed(3, matrix(stats::rnorm(60 * 3, 0, .01), 60, 3))

  for (lam in c(0.2, 0.02)) {
    ref <- .fit_sgl_r(X, Y, lam, .5, B_init = B0, max_iter = 300, tol = 1e-6)
    expect_equal(.fit_sgl(X, Y, lam, .5, B_init = B0, max_iter = 300, tol = 1e-6),
                 ref, tolerance = 1e-10)
    expect_equal(.cd_solve(cd_gram, 0L, .sgl_pars(lam, .5), B0 = B0,
                           max_iter = 300, tol = 1e-6),
                 ref, tolerance = 1e-10)
  }
})

test_that("multivariate_sgl gives the same fit with and without X'X held", {
  X <- sim_geno(80, 30, 4)
  Y <- sim_expr(X, 3, 5)
  a <- multivariate_sgl(X, Y, nlambda = 6, nfolds = 3, seed = 9)
  b <- withr::with_options(list(isotwas.gram_max_p = 0),
                           multivariate_sgl(X, Y, nlambda = 6, nfolds = 3, seed = 9))
  expect_equal(a$summary$r2, b$summary$r2, tolerance = 1e-10)
  expect_equal(get_weight_matrix(a), get_weight_matrix(b), tolerance = 1e-10)
})
