# Tests for the scoring and selection helpers in R/utils.R.
#
# Several tests here are CHARACTERIZATION tests: they record what the package
# does today, including behaviour that is arguably wrong, so that any change is
# a deliberate, visible edit to this file rather than a silent drift. Each one
# says so and explains what the alternative would be.

# ---- calc.mse ----------------------------------------------------------

test_that("calc.mse returns one mean squared error per column", {

    obs <- cbind(a = c(1, 2, 3), b = c(1, 1, 1))
    pred <- cbind(a = c(1, 2, 3), b = c(2, 2, 2))

    expect_equal(unname(calc.mse(obs, pred)), c(0, 1))
})

test_that("calc.mse accepts vectors", {
    expect_equal(unname(calc.mse(c(1, 2, 3), c(1, 2, 4))), 1 / 3)
})

# ---- calc.r2 -----------------------------------------------------------

test_that("calc.r2 returns R2 and P for the i-th column", {

    set.seed(11)
    obs <- cbind(stats::rnorm(50))
    pred <- cbind(obs[, 1] + stats::rnorm(50, sd = 0.1))

    got <- calc.r2(1, obs, pred)

    expect_named(got, c("R2", "P"))
    expect_gt(got$R2, 0.9)
    expect_lt(got$P, 1e-10)
})

test_that("calc.r2 scores an uninformative prediction near zero or below", {

    set.seed(12)
    obs <- cbind(stats::rnorm(50))
    pred <- cbind(stats::rnorm(50))

    expect_lt(calc.r2(1, obs, pred)$R2, 0.1)
})

test_that("calc.r2 maps a degenerate fit to 0 rather than NA", {

    # A constant prediction gives an NA adjusted R2 from lm(); the function
    # maps that to 0 explicitly.
    obs <- cbind(c(1, 2, 3, 4))
    pred <- cbind(c(2, 2, 2, 2))

    expect_warning(got <- calc.r2(1, obs, pred))
    expect_identical(got$R2, 0)
})

test_that("calc.r2 refits, so it is invariant to rescaling the prediction", {

    # This is the property that makes calc.r2 incomparable with a direct
    # out-of-sample R2: it regresses obs on pred, so any affine transform of
    # pred scores identically. A method that scores itself with
    # 1 - ss_res/ss_tot is penalised for miscalibration where this is not.
    set.seed(13)
    obs <- cbind(stats::rnorm(60))
    pred <- cbind(obs[, 1] * 0.5 + stats::rnorm(60, sd = 0.4))

    plain <- calc.r2(1, obs, pred)$R2
    rescaled <- calc.r2(1, obs, cbind(pred[, 1] * 3 + 10))$R2

    expect_equal(plain, rescaled)
})

# ---- get_best_isotwas --------------------------------------------------

test_that("selection is per transcript, not per gene", {

    # t1 is best under method A, t2 under method B. The assembled model must
    # mix them rather than picking one method for the whole gene.
    a <- mk_model("A", c(0.50, 0.10))
    b <- mk_model("B", c(0.20, 0.80))

    got <- get_best_isotwas(list(a, b), G = 2)

    expect_identical(got$best_methods$best_method, c("A", "B"))
    expect_equal(got$transcripts[["t1"]]$r2, 0.50)
    expect_equal(got$transcripts[["t2"]]$r2, 0.80)
})

test_that("the assembled model is labelled 'best' and carries the R2 matrix", {

    a <- mk_model("A", c(0.5, 0.1))
    b <- mk_model("B", c(0.2, 0.8))

    got <- get_best_isotwas(list(a, b), G = 2)

    expect_s3_class(got, "isotwas_model")
    expect_identical(got$method, "best")
    expect_identical(dim(got$r2_comparison), c(2L, 2L))
    expect_identical(colnames(got$r2_comparison), c("A", "B"))
    expect_identical(rownames(got$r2_comparison), c("t1", "t2"))
})

test_that("a transcript missing from a model scores empty_r2 for that model", {

    a <- mk_model("A", c(0.5, 0.1))
    b <- mk_model("B", c(0.2, 0.8))
    b$transcripts[["t2"]] <- NULL

    got <- get_best_isotwas(list(a, b), G = 2)

    expect_equal(got$r2_comparison["t2", "B"], -Inf)
    expect_identical(got$best_methods$best_method[2], "A")

    legacy <- get_best_isotwas(list(a, b), G = 2, empty_r2 = 0)
    expect_equal(legacy$r2_comparison["t2", "B"], 0)
})

test_that("an NA R2 is scored as empty_r2", {

    a <- mk_model("A", c(NA_real_, 0.1))
    b <- mk_model("B", c(0.2, 0.8))

    got <- get_best_isotwas(list(a, b), G = 2)
    expect_equal(got$r2_comparison["t1", "A"], -Inf)

    legacy <- get_best_isotwas(list(a, b), G = 2, empty_r2 = 0)
    expect_equal(legacy$r2_comparison["t1", "A"], 0)
})

test_that("empty_r2 = 0 keeps the historical scoring of an empty model", {

    # An empty model scored 0 outranks a real model whose adjusted R2 is
    # negative, and is then emitted as that transcript's weights. This is what
    # r2_scale = 'legacy' uses.
    empty <- mk_model("A", c(0.40), n_snps = 0L)
    real <- mk_model("B", c(-0.01), n_snps = 30L)

    got <- get_best_isotwas(list(empty, real), G = 1, empty_r2 = 0)

    expect_equal(got$r2_comparison["t1", "A"], 0)
    expect_identical(got$best_methods$best_method, "A")
    expect_identical(nrow(got$transcripts[["t1"]]$weights), 0L)
})

test_that("by default a real model beats an empty one", {

    empty <- mk_model("A", c(0.40), n_snps = 0L)
    real <- mk_model("B", c(-0.01), n_snps = 30L)

    got <- get_best_isotwas(list(empty, real), G = 1)

    expect_identical(got$best_methods$best_method, "B")
    expect_identical(nrow(got$transcripts[["t1"]]$weights), 30L)
})

test_that("a transcript no model could fit is labelled NA, not given a winner", {

    # With every candidate empty, which.max would otherwise return the first
    # model and emit its empty weights as though they had been selected.
    a <- mk_model("A", c(0.40), n_snps = 0L)
    b <- mk_model("B", c(0.20), n_snps = 0L)

    got <- get_best_isotwas(list(a, b), G = 1, empty_r2 = -Inf)

    expect_true(is.na(got$best_methods$best_method))
    expect_true(is.na(got$transcripts[["t1"]]$method))
})

test_that("a nested ensemble's leaf label survives selection", {

    # compute_isotwas builds its `univariate` entry by calling get_best() a
    # second time, so that entry's own $method is the literal "best". Stamping
    # the parent's $method destroyed which of enet/BLUP/SuSiE won -- making
    # the reported best method a constant for every transcript the ensemble
    # took.
    flat <- mk_model("spls", c(0.10, 0.90))
    ens <- mk_ensemble(c(0.50, 0.40),
                       leaves = c("univariate_blup", "univariate_susie"))

    got <- get_best_isotwas(list(flat, ens), G = 2)

    expect_identical(got$best_methods$best_method,
                     c("univariate_blup", "spls"))
    expect_identical(got$transcripts[["t1"]]$method, "univariate_blup")
    expect_false(any(got$best_methods$best_method == "best"))
})

test_that("resolving the label does not change which model wins", {

    flat <- mk_model("spls", c(0.10, 0.90))
    ens <- mk_ensemble(c(0.50, 0.40),
                       leaves = c("univariate_blup", "univariate_susie"))

    got <- get_best_isotwas(list(flat, ens), G = 2)

    # t1 from the ensemble, t2 from spls -- the weights, not just the names.
    expect_equal(got$transcripts[["t1"]]$r2, 0.50)
    expect_equal(got$transcripts[["t2"]]$r2, 0.90)
    expect_identical(colnames(got$r2_comparison), c("spls", "best"))
})

test_that(".leaf_method falls back for a flat model or unknown transcript", {

    flat <- mk_model("spls", c(0.1, 0.2))
    expect_identical(.leaf_method(flat, "t1", "spls"), "spls")

    ens <- mk_ensemble(c(0.5, 0.4), leaves = c("univariate_blup", "univariate_susie"))
    expect_identical(.leaf_method(ens, "t1", "best"), "univariate_blup")
    expect_identical(.leaf_method(ens, "nope", "best"), "best")
})

test_that(".score_r2 reproduces both scales in use", {

    set.seed(61)
    obs <- stats::rnorm(60)
    pred <- obs * 0.5 + stats::rnorm(60, sd = 0.4)
    n <- length(obs)

    # "refit" is calc.r2()
    expect_equal(.score_r2(obs, pred, "refit")$R2,
                 summary(stats::lm(obs ~ pred))$adj.r.sq)

    # "oos" is the formula compute_mrce() and multivariate_sgl() inline
    ss_res <- sum((obs - pred)^2)
    ss_tot <- sum((obs - mean(obs))^2)
    expect_equal(.score_r2(obs, pred, "oos")$R2,
                 1 - (1 - (1 - ss_res / ss_tot)) * (n - 1) / (n - 2))
})

test_that("the two R2 scales are not interchangeable", {

    # This is the incompatibility, stated as an executable fact: "refit"
    # regresses obs on pred, so it cannot see a miscalibrated prediction,
    # while "oos" is penalised for one. Mixing them in a single which.max
    # compares different quantities.
    set.seed(62)
    obs <- stats::rnorm(60)
    pred <- obs * 0.5 + stats::rnorm(60, sd = 0.4)

    expect_equal(.score_r2(obs, pred, "refit")$R2,
                 .score_r2(obs, pred * 3 + 10, "refit")$R2)
    expect_false(isTRUE(all.equal(.score_r2(obs, pred, "oos")$R2,
                                  .score_r2(obs, pred * 3 + 10, "oos")$R2)))
})

test_that("an unpredictable transcript scores negative under 'oos', not 0", {

    # Flooring this at 0 is what made a floored method outrank every
    # calc.r2()-scored rival on transcripts nothing can predict.
    set.seed(63)
    obs <- stats::rnorm(80)
    noise <- stats::rnorm(80)

    expect_lt(.score_r2(obs, noise, "oos")$R2, 0)
    expect_equal(max(0, .score_r2(obs, noise, "oos")$R2), 0)
})

test_that("flooring one method's score biases selection toward it", {

    # End to end over a transcript no method can predict: the floored method
    # reports 0, the calc.r2-scored method reports a negative number, and
    # which.max takes the floored one despite it selecting nothing useful.
    floored <- mk_model("mrce_lasso", c(0), n_snps = 5L)
    honest <- mk_model("spls", c(-0.02), n_snps = 12L)

    got <- get_best_isotwas(list(floored, honest), G = 1, empty_r2 = 0)
    expect_identical(got$best_methods$best_method, "mrce_lasso")

    # With both scored on the same scale the comparison is between two
    # negatives and the better one wins.
    floored_unclamped <- mk_model("mrce_lasso", c(-0.15), n_snps = 5L)
    got2 <- get_best_isotwas(list(floored_unclamped, honest), G = 1)
    expect_identical(got2$best_methods$best_method, "spls")
})
