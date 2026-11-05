# .susie_coef() must reproduce susieR::susie() as univariate_susie called it.

test_that(".susie_coef matches susieR::susie as univariate_susie calls it", {
  for (dims in list(c(150, 50), c(50, 150))) {
    X <- sim_geno(dims[1], dims[2], 16)
    X[, 2:3] <- X[, 1]
    Y <- sim_expr(X, 2, 17) + 3
    ref <- sapply(1:2, function(g) {
      fit <- suppressMessages(suppressWarnings(susieR::susie(
        Matrix::Matrix(X, sparse = TRUE), Y[, g],
        estimate_residual_variance = TRUE, estimate_prior_variance = FALSE,
        scaled_prior_variance = 0.1, verbose = FALSE)))
      as.numeric(stats::coef(fit))[-1]
    })
    new <- suppressWarnings(.susie_coef(X, Y))
    expect_equal(new, ref, tolerance = 1e-8, ignore_attr = TRUE)
  }
})
