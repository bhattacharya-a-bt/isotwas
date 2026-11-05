# multivariate_mtlasso's compiled solver must reproduce the R solver it
# replaces.

test_that("the compiled L21 multi-task lasso matches .fit_mtl_l21_r", {
  X <- scale(sim_geno(120, 60, 1), scale = FALSE)
  Y <- scale(sim_expr(X, 3, 2), scale = FALSE)
  cd_gram <- .cd_data(X, Y, NULL, crossprod(X), crossprod(X, Y))
  B0 <- withr::with_seed(3, matrix(stats::rnorm(60 * 3, 0, .01), 60, 3))

  for (lam in c(0.2, 0.02)) {
    ref <- .fit_mtl_l21_r(X, Y, lam, W_init = B0, max_iter = 300, tol = 1e-6)
    expect_equal(.fit_mtl_l21(X, Y, lam, W_init = B0, max_iter = 300, tol = 1e-6),
                 ref, tolerance = 1e-10)
    expect_equal(.cd_solve(cd_gram, 1L, c(lam, 0, 0), B0 = B0,
                           max_iter = 300, tol = 1e-6),
                 ref, tolerance = 1e-10)
  }
})
