# compute_isotwas compares methods by a which.max of cross-validated R2, so
# every method must report the same definition of it. r2_scale = 'coherent'
# (the default) scores every method by the adjusted R2 of lm(obs ~ pred);
# 'legacy' keeps the out-of-sample R2 of mrce_lasso, sgl, mtlasso, stacking
# and graph_reg, clamped at 0.

test_that("on the same prediction 'oos' is never above 'refit'", {
  set.seed(261)
  obs <- stats::rnorm(80)
  noise <- stats::rnorm(80)
  preds <- list(good = obs + noise,
                shrunk = 0.2 * (obs + noise),
                shifted = obs + noise + 3,
                stretched = 4 * (obs + noise) - 1,
                weak = 0.1 * obs + noise,
                flipped = -(obs + noise),
                constant = rep(0.5, 80))
  for (p in preds) {
    expect_lte(.score_r2(obs, p, "oos")$R2, .score_r2(obs, p, "refit")$R2)
  }
  # "refit" is unchanged by the scale and offset of a prediction
  expect_equal(.score_r2(obs, preds$shrunk)$R2, .score_r2(obs, preds$good)$R2)
  expect_equal(.score_r2(obs, preds$shifted)$R2, .score_r2(obs, preds$good)$R2)
})

test_that("a mis-scaled but well-correlated prediction no longer loses", {

  # Method A's prediction tracks expression closely but is shrunk and
  # shifted, as a penalised fit without its intercept is; method B's is
  # well scaled but noisier. Scored by "oos", A loses to B scored by
  # "refit", although A is the better predictor of what a TWAS tests.
  set.seed(262)
  obs <- stats::rnorm(100)
  pred_a <- 0.3 * (obs + stats::rnorm(100, sd = 0.5)) + 2
  pred_b <- obs + stats::rnorm(100, sd = 1.2)
  expect_gt(stats::cor(obs, pred_a), stats::cor(obs, pred_b))

  pick <- function(r2_a, r2_b) {
    a <- mk_model("mrce_lasso", r2_a, n_snps = 5L)
    b <- mk_model("multi_enet", r2_b, n_snps = 5L)
    get_best_isotwas(list(a, b), G = 1)$best_methods$best_method
  }
  refit_b <- .score_r2(obs, pred_b, "refit")$R2

  # legacy: A scored by "oos" and clamped, B by "refit"
  expect_identical(pick(max(0, .score_r2(obs, pred_a, "oos")$R2), refit_b),
                   "multi_enet")
  # coherent: both scored by "refit", the default of every fitter
  expect_identical(pick(.score_r2(obs, pred_a)$R2, refit_b), "mrce_lasso")
})

test_that("compute_isotwas gives every method the R2 of its r2_scale", {
  skip_if_not_installed("huge")
  X <- sim_geno(150, 30, 263)
  Y <- sim_expr(X, 3, 264)
  S <- matrix(.3, 3, 3); diag(S) <- 1
  run <- function(scale) {
    compute_isotwas(X, Y, Y.rep = Y, R = 1, id = seq_len(nrow(Y)),
                    omega_est = "mean", omega_nlambda = 3,
                    predict_nlambda = 6, sgl_nlambda = 6,
                    similarity_matrix = S, r2_scale = scale,
                    verbose = FALSE, seed = 265, return_all = TRUE)
  }
  oos_methods <- c("mrce_lasso", "sgl", "mtlasso", "stacking", "graph_reg")
  for (scale in c("coherent", "legacy")) {
    res <- run(scale)
    # mrce_lasso (the best fit over Omegas) and the univariate ensemble are
    # both labelled "best"; they come second and last
    labels <- vapply(res$method_results, function(m) m$method, "")
    expect_identical(which(labels == "best"), c(2L, length(labels)))
    labels[labels == "best"] <- c("mrce_lasso", "univariate")
    expect_setequal(labels, isotwas_methods())
    for (i in seq_along(res$method_results)) {
      m <- res$method_results[[i]]
      label <- labels[i]
      if (scale == "legacy" && label %in% oos_methods) {
        expect_scored_by(m, Y, "oos", floor = TRUE)
      } else {
        expect_scored_by(m, Y, "refit")
      }
      # the comparison table reports what each method scored
      expect_equal(res$comparison[[label]], unname(m$summary$r2))
    }
  }
})
