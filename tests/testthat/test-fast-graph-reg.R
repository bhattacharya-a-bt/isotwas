# multivariate_graph_reg's compiled solver must reproduce the R solver it
# replaces.

test_that("the compiled graph-regularised solver matches .fit_graph_reg_r", {
  X <- scale(sim_geno(120, 60, 1), scale = FALSE)
  Y <- scale(sim_expr(X, 3, 2), scale = FALSE)
  B0 <- withr::with_seed(3, matrix(stats::rnorm(60 * 3, 0, .01), 60, 3))
  S <- matrix(.3, 3, 3); diag(S) <- 1
  L <- .compute_laplacian(S)
  for (lam in c(0.2, 0.02)) {
    for (lg in c(0, 0.4)) {
      ref <- .fit_graph_reg_r(X, Y, L, lam, lg, .5, B_init = B0,
                              max_iter = 300, tol = 1e-6)
      expect_equal(.fit_graph_reg(X, Y, L, lam, lg, .5, B_init = B0,
                                  max_iter = 300, tol = 1e-6),
                   ref, tolerance = 1e-10)
    }
  }
})
