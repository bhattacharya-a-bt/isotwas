# The compiled coordinate-descent solver gets each row's gradient either from
# X'X (covariance mode) or from the residuals (naive mode, used when X'X is
# too large to hold). The two make the same updates, so they must agree.

test_that("covariance and naive modes of the compiled solver agree", {
  X <- scale(sim_geno(120, 60, 1), scale = FALSE)
  Y <- scale(sim_expr(X, 3, 2), scale = FALSE)
  idx <- 1:90
  gram <- .cd_data(X, Y, idx, crossprod(X), crossprod(X, Y))
  naive <- .cd_data(X, Y, idx)
  expect_null(naive$G)
  # the training set's X'X by subtraction equals a fresh crossprod
  expect_equal(gram$G, crossprod(X[idx, ]), tolerance = 1e-12)

  B0 <- withr::with_seed(3, matrix(stats::rnorm(60 * 3, 0, .01), 60, 3))
  S <- matrix(.3, 3, 3); diag(S) <- 1
  L <- diag(rowSums(S)) - S
  e <- eigen(L, symmetric = TRUE)
  for (lam in c(0.2, 0.02)) {
    for (type in 0:3) {
      pars <- c(lam * .5, lam * .5, if (type == 3) .4 else 0)
      a <- .cd_solve(gram, type, pars, B0 = B0, max_iter = 300, tol = 1e-8,
                     V = e$vectors, ev = e$values)
      b <- .cd_solve(naive, type, pars, B0 = B0, max_iter = 300, tol = 1e-8,
                     V = e$vectors, ev = e$values)
      expect_equal(a, b, tolerance = 1e-10)
      expect_gt(sum(a != 0), 0)
    }
  }
})
