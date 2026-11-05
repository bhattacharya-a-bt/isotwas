# blasso() never updates an entry's own soft-threshold value, so rblasso()
# must start each one without that entry's own term. Otherwise a warm start
# carries an offset for the whole fit and converges somewhere else.

test_that("a warm-started MRCE lasso converges to the cold-start solution", {
  withr::local_seed(6)
  X <- scale(matrix(stats::rbinom(150 * 40, 2, .3), 150, 40), scale = FALSE)
  Y <- scale(X[, 1:4] * .4 + matrix(stats::rnorm(150 * 4), 150, 4),
             scale = FALSE)
  Om <- solve(stats::cov(Y)); Om <- (Om + t(Om)) / 2
  S <- crossprod(X); M <- crossprod(X, Y) %*% Om
  # the stopping rule is on the objective, quadratic in the coefficient
  # change, so coefficients settle to about sqrt(tol)
  tol <- 1e-20 * sum(crossprod(Y) * Om)
  nl <- matrix(150 * .02, 40, 4)
  B0 <- matrix(0, 40, 4)
  B0[sample(160, 20)] <- .3

  # the compiled solver (default) and the reference blasso()
  for (ref in c(FALSE, TRUE)) {
    withr::local_options(list(isotwas.mrce_reference = ref))
    cold <- rblasso(S, M, Om, nl, 150, tol = tol, maxit = 5000)
    expect_gt(sum(cold$B != 0), 0)
    # restarting at the solution must not move it
    again <- rblasso(S, M, Om, nl, 150, B0 = cold$B, tol = tol, maxit = 5000)
    expect_equal(again$B, cold$B, tolerance = 1e-8)
    # nor may an arbitrary starting point change where it ends up
    warm <- rblasso(S, M, Om, nl, 150, B0 = B0, tol = tol, maxit = 5000)
    expect_equal(warm$B, cold$B, tolerance = 1e-8)
  }
})
