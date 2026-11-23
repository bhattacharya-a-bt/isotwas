# univariate_blup is ridge BLUP: rrBLUP::mixed.solve(y, Z = X, K = diag(m)).
# These tests check the fast solver against mixed.solve and the behaviour of
# the weights.

# Weights are "dense same-sign" when at least 95% share one sign and
# |sum(w)| / sum(|w|) is at least 0.95.
same_sign <- function(w) {
  max(mean(w > 0), mean(w < 0)) >= .95 && abs(sum(w)) / sum(abs(w)) >= .95
}

test_that(".blup_fit matches rrBLUP::mixed.solve with K = I in both spectral branches", {
  skip_if_not_installed("rrBLUP")
  # (120, 40): mixed.solve's "cholesky" branch, .blup_prep's X'X;
  # (41, 40): mixed.solve's "eigen" branch, .blup_prep's X'X;
  # (40, 120): "eigen" and XX'
  for (dims in list(c(120, 40), c(41, 40), c(40, 120))) {
    X <- sim_geno(dims[1], dims[2], 14)
    X[, 2:4] <- X[, 1]   # duplicated SNPs make the design rank deficient
    X[, 5] <- 1          # and a monomorphic one
    Y <- sim_expr(X, 2, 15) + 5
    Xc <- sweep(X, 2, colMeans(X))
    new <- .blup_fit(Y, X)
    new_fit <- sweep(X %*% new$weights, 2, new$intercept, "+")

    # on centred dosages, and on the raw dosages, where mixed.solve's beta is
    # the intercept for the raw-scale weights
    for (Z in list(Xc, X)) {
      ref <- lapply(1:2, function(g)
        rrBLUP::mixed.solve(Y[, g], Z = Z, K = diag(ncol(X))))
      ref_u <- sapply(ref, function(r) as.numeric(r$u))
      ref_fit <- sapply(ref, function(r) drop(r$beta) + drop(Z %*% r$u))
      # mixed.solve finds the REML optimum with optimize()'s default
      # tolerance, so the weights agree to about 1e-4 relative at worst
      expect_lt(max(abs(new$weights - ref_u)) / max(abs(ref_u)), 1e-4)
      expect_lt(max(abs(new_fit - ref_fit)) /
                  max(abs(sweep(ref_fit, 2, colMeans(ref_fit)))), 1e-5)
    }
    expect_equal(unname(new$intercept), sapply(ref, function(r) drop(r$beta)),
                 tolerance = 1e-4)
    expect_identical(new$weights[5, ], c(0, 0))
  }
})

test_that(".blup_fit fits responses with missing values on their observed rows", {
  X <- sim_geno(60, 30, 16)
  Y <- sim_expr(X, 4, 17)
  Y[c(3, 10), 2:3] <- NA
  Y[, 4] <- NA
  fit <- .blup_fit(Y, X)
  ok <- !is.na(Y[, 2])
  expect_equal(fit$weights[, 1], .blup_fit(Y[, 1], X)$weights[, 1])
  expect_equal(fit$weights[, 2:3], .blup_fit(Y[ok, 2:3], X[ok, ])$weights)
  # no observed values: no fit
  expect_identical(fit$weights[, 4], rep(0, 30))
  expect_true(is.na(fit$intercept[4]))
  expect_error(.blup_fit(Y[, 1], replace(X, 1, NA)), "missing values")
  expect_error(.blup_fit(Y[0, 1], X[0, ]), "no samples")
})

test_that("BLUP weights are not one-signed under the null", {
  # Under the old prior (K = X'X / (m - 1) on uncentred dosages) about 60% of
  # fits like these were dense same-sign, with |cor(w, ALT frequency)|
  # about 0.8.
  res <- withr::with_seed(42, t(replicate(20, {
    p <- stats::runif(400, .05, .5)
    X <- sapply(p, function(f) stats::rbinom(208, 2, f))
    y <- 7 + stats::rnorm(208)
    w <- drop(.blup_fit(y, X)$weights)
    c(same_sign = same_sign(w), cor_p = stats::cor(w, p),
      pos = mean(w > 0))
  })))
  expect_equal(sum(res[, "same_sign"]), 0)
  expect_lt(stats::median(abs(res[, "cor_p"])), .2)
  expect_true(all(res[, "pos"] > .3 & res[, "pos"] < .7))
})

test_that("recoding a SNP's alleles flips only that SNP's weight", {
  X <- sim_geno(150, 60, 18)
  Y <- sim_expr(X, 1, 19)
  flip <- c(2, 7, 30)
  X2 <- X
  X2[, flip] <- 2 - X2[, flip]
  w <- drop(.blup_fit(Y, X)$weights)
  w2 <- drop(.blup_fit(Y, X2)$weights)
  # exact in exact arithmetic; each fit reruns the REML search
  expect_equal(w2[flip], -w[flip], tolerance = 1e-6)
  expect_equal(w2[-flip], w[-flip], tolerance = 1e-6)
})

test_that("BLUP predicts simulated cis signal out of sample", {
  # 5 causal SNPs of random sign, cis-h2 = 0.1, n = 200 training and 1000
  # test samples: correlation of the prediction with the true genetic value.
  # Ridge BLUP averages about 0.26 here (the old prior about 0.09, with some
  # replicates below 0).
  r <- withr::with_seed(7, replicate(10, {
    n <- 200
    m <- 300
    p <- stats::runif(m, .05, .5)
    Xa <- sapply(p, function(f) stats::rbinom(n + 1000, 2, f))
    beta <- numeric(m)
    beta[sample(m, 5)] <- stats::rnorm(5)
    g <- drop(Xa %*% beta)
    g <- g / stats::sd(g[1:n]) * sqrt(.1)
    y <- g[1:n] + stats::rnorm(n, sd = sqrt(.9)) + 7
    w <- .blup_fit(y, Xa[1:n, ])$weights
    stats::cor(drop(Xa[-(1:n), ] %*% w), g[-(1:n)])
  }))
  expect_gt(mean(r), .2)
  expect_true(all(r > 0))
})

test_that("univariate_blup's cross-validated predictions carry the fold intercepts", {
  X <- sim_geno(100, 50, 20)
  Y <- sim_expr(X, 2, 21)
  fit <- function(Y) univariate_blup(X, Y, diag(2), nfolds = 4,
                                     verbose = FALSE, seed = 3)
  a <- fit(Y)
  b <- fit(Y + 100)
  for (tx in names(a$transcripts)) {
    ta <- a$transcripts[[tx]]
    tb <- b$transcripts[[tx]]
    # adding a constant moves the predictions, not the weights or R2
    expect_equal(tb$weights, ta$weights, tolerance = 1e-6)
    expect_equal(tb$predicted, ta$predicted + 100, tolerance = 1e-6)
    expect_equal(tb$r2, ta$r2, tolerance = 1e-6)
  }

  # each held-out prediction is X.test %*% w + intercept from its fold's fit
  folds <- withr::with_seed(3, caret::createFolds(1:100, k = 4, returnTrain = TRUE))
  pred <- matrix(NA_real_, 100, 2)
  for (tr in folds) {
    f <- .blup_fit(Y[tr, ], X[tr, ])
    pred[-tr, ] <- sweep(X[-tr, ] %*% f$weights, 2, f$intercept, "+")
  }
  expect_equal(unname(sapply(a$transcripts, `[[`, "predicted")), pred)
  expect_equal(unname(sapply(a$transcripts, `[[`, "r2")),
               sapply(1:2, function(g) .score_r2(Y[, g], pred[, g])$R2))
})

test_that("compute_blup is deprecated and returns univariate_blup's fit", {
  X <- sim_geno(80, 30, 22)
  Y <- sim_expr(X, 2, 23)
  expect_warning(old <- compute_blup(X, Y, diag(2), nfolds = 4,
                                     verbose = FALSE, seed = 3),
                 "deprecated")
  new <- univariate_blup(X, Y, diag(2), nfolds = 4, verbose = FALSE, seed = 3)
  for (i in 1:2) {
    tm <- new$transcripts[[i]]
    expect_identical(old[[i]]$Transcript, tm$transcript_id)
    expect_identical(old[[i]]$Model, tm$weights)
    expect_identical(old[[i]]$Pred, tm$predicted)
    expect_identical(old[[i]]$R2, tm$r2)
  }
})
