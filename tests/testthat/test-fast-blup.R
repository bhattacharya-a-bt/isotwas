# .blup_u() must reproduce rrBLUP::mixed.solve() as univariate_blup called it.

test_that(".blup_u matches rrBLUP::mixed.solve in both spectral branches", {
  for (dims in list(c(120, 40), c(40, 120))) {
    X <- sim_geno(dims[1], dims[2], 14)
    X[, 2:4] <- X[, 1]   # duplicated SNPs make the design rank deficient
    Y <- sim_expr(X, 2, 15)
    K <- crossprod(X) / (ncol(X) - 1)
    ref <- sapply(1:2, function(g)
      as.numeric(rrBLUP::mixed.solve(Y[, g], Z = X, K = K)$u))
    new <- .blup_u(Y, X)
    # mixed.solve finds the REML optimum with optimize()'s default tolerance,
    # about 1e-4 relative, so the two agree to that precision
    expect_lt(max(abs(new - ref)) / max(abs(ref)), 1e-4)
  }
})
