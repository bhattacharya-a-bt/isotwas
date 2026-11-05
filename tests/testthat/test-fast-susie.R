# .susie_coef() must reproduce susieR 0.14.2's susie() as univariate_susie
# called it, whatever susieR version is installed.

test_that(".susie_coef matches susieR::susie 0.14.x as univariate_susie calls it", {
  skip_if_not(utils::packageVersion("susieR") >= "0.14.0" &&
                utils::packageVersion("susieR") < "0.15.0",
              "susieR 0.14.x is not installed; the fixed values below still apply")
  for (dims in list(c(150, 50), c(50, 150))) {
    X <- sim_geno(dims[1], dims[2], 16)
    X[, 2:3] <- X[, 1]
    X[, 4] <- 2  # no variance
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

test_that(".susie_coef gives susieR 0.14.2's fit with any susieR installed", {
  # coef(susieR::susie(...))[-1] under susieR 0.14.2, as above
  expected <- list(
    c(0.001436254447, 0.001436254447, 0, -0.001877321659, 0.9870244069,
      -0.004655664581, -0.08341767274, 0.001217455013, -0.02202901031,
      0.07378357395),
    c(0.01102747326, 0.01102747326, 0, 0.07398275954, -0.08967264122,
      0.06918883785, -0.05481534304, -0.03243388214, -0.02575505246,
      0.06334463516, 0.2381715737, -0.005150114802))
  dims <- list(c(60, 10), c(8, 12))  # n > p uses X'X, n < p fitted values
  for (k in 1:2) {
    X <- sim_geno(dims[[k]][1], dims[[k]][2], 21)
    X[, 2] <- X[, 1]
    X[, 3] <- 1
    y <- sim_expr(X, 1, 22, n_causal = 2)[, 1] + 3
    new <- suppressWarnings(.susie_coef(X, cbind(y)))[, 1]
    expect_equal(new, expected[[k]], tolerance = 1e-8, ignore_attr = TRUE)
  }
})
