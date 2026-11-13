# compute_mrce fits on centred data, so its cross-validated predictions need
# each fold's intercept, my - mx'B. Without it they were off by roughly the
# expression mean, which "oos" scoring penalised in full.

test_that("compute_mrce's CV predictions carry the fold intercept", {
  X <- sim_geno(150, 30, 281)
  Y <- sim_expr(X, 3, 282)
  Om <- solve(stats::cov(Y))
  fit <- function(Y, ...) compute_mrce(X, Y, nlambda = 8, Omega = Om,
                                       nfolds = 5, tol.in = 1e-3,
                                       maxit.in = 100, verbose = FALSE,
                                       seed = 283, ...)
  shift <- 5
  for (m in c("refit", "oos")) {
    a <- fit(Y, r2_method = m)
    b <- fit(Y + shift, r2_method = m)
    # Adding a constant to expression shifts every prediction by it and
    # changes nothing else: the lambda path, the choice, the weights, the R2
    expect_equal(a$lambda, b$lambda, tolerance = 1e-12)
    expect_identical(match(a$best_lambda, a$lambda),
                     match(b$best_lambda, b$lambda))
    expect_equal(a$cv_r2, b$cv_r2, tolerance = 1e-8)
    for (tx in names(a$transcripts)) {
      expect_equal(b$transcripts[[tx]]$predicted,
                   a$transcripts[[tx]]$predicted + shift, tolerance = 1e-8)
      expect_equal(b$transcripts[[tx]]$weights, a$transcripts[[tx]]$weights,
                   tolerance = 1e-8)
      expect_equal(b$transcripts[[tx]]$r2, a$transcripts[[tx]]$r2,
                   tolerance = 1e-8)
    }
  }
  # Held-out predictions are centred near the expression mean
  a <- fit(Y + shift)
  for (tx in names(a$transcripts)) {
    expect_lt(abs(mean(a$transcripts[[tx]]$predicted) - mean(Y[, tx] + shift)),
              0.1 * stats::sd(Y[, tx]))
  }
})
