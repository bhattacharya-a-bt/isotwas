# The compiled MRCE kernel must reproduce blasso(), and compute_mrce() must give
# the same fit with a fold cache shared across Omegas.

test_that("the fast MRCE solver reproduces blasso, warm starts included", {
  X <- scale(sim_geno(150, 40, 6), scale = FALSE)
  Y <- scale(sim_expr(X, 4, 7), scale = FALSE)
  Om <- solve(stats::cov(Y)); Om <- (Om + t(Om)) / 2
  S <- crossprod(X); M <- crossprod(X, Y) %*% Om
  tol <- 1e-9 * sum(crossprod(Y) * Om)
  B0 <- matrix(0, 40, 4)
  B0[withr::with_seed(8, sample(160, 20))] <- .05
  for (lam in c(.2, .01)) {
    nl <- matrix(150 * lam, 40, 4)
    for (b0 in list(NULL, B0)) {
      ref <- withr::with_options(list(isotwas.mrce_reference = TRUE),
                                 rblasso(S, M, Om, nl, 150, B0 = b0, tol = tol, maxit = 500))
      new <- rblasso(S, M, Om, nl, 150, B0 = b0, tol = tol, maxit = 500)
      expect_equal(new$B, ref$B, tolerance = 1e-10)
      expect_equal(new$iterations, ref$iterations)
    }
  }
})

test_that("compute_mrce gives the same fit with a shared fold cache", {
  X <- sim_geno(80, 25, 10)
  Y <- sim_expr(X, 3, 11)
  Om <- solve(stats::cov(Y))
  a <- compute_mrce(X, Y, nlambda = 5, Omega = Om, nfolds = 3, tol.in = 1e-3,
                    maxit.in = 100, verbose = FALSE, seed = 12)
  b <- compute_mrce(X, Y, nlambda = 5, Omega = Om, nfolds = 3, tol.in = 1e-3,
                    maxit.in = 100, verbose = FALSE, seed = 12,
                    fold_cache = .mrce_fold_cache(X, Y, 3, 12))
  expect_equal(a$summary$r2, b$summary$r2)
  expect_equal(get_weight_matrix(a), get_weight_matrix(b))
  expect_error(compute_mrce(X, Y, nlambda = 5, Omega = Om, nfolds = 3,
                            tol.in = 1e-3, maxit.in = 100, verbose = FALSE,
                            seed = 13, fold_cache = .mrce_fold_cache(X, Y, 3, 12)),
               "different folds")
})
