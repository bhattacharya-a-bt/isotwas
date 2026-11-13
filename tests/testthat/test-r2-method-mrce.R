# compute_mrce is the one method whose fit, not only its reported score,
# follows r2_method: it chooses lambda by the mean of that score across
# transcripts.

test_that("compute_mrce chooses lambda by the mean of its r2_method score", {
  X <- sim_geno(150, 30, 307)
  Y <- sim_expr(X, 3, 308)
  Om <- solve(stats::cov(Y))
  fit <- function(...) compute_mrce(X, Y, nlambda = 8, Omega = Om, nfolds = 5,
                                    tol.in = 1e-3, maxit.in = 100,
                                    verbose = FALSE, seed = 309, ...)
  a <- fit()
  b <- fit(r2_method = "oos")
  expect_identical(a$lambda, b$lambda)
  for (m in list(a, b)) {
    expect_identical(m$best_lambda, m$lambda[which.max(colMeans(m$cv_r2))])
  }
  # Here "oos", which also penalises the scale of a prediction, picks a
  # larger lambda than "refit"
  expect_gt(b$best_lambda, a$best_lambda)
})
