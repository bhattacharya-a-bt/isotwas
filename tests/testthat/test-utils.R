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

test_that("a transcript missing from a model scores 0 for that model", {

    a <- mk_model("A", c(0.5, 0.1))
    b <- mk_model("B", c(0.2, 0.8))
    b$transcripts[["t2"]] <- NULL

    got <- get_best_isotwas(list(a, b), G = 2)

    expect_equal(got$r2_comparison["t2", "B"], 0)
    expect_identical(got$best_methods$best_method[2], "A")
})

test_that("an NA R2 is treated as 0", {

    a <- mk_model("A", c(NA_real_, 0.1))
    b <- mk_model("B", c(0.2, 0.8))

    got <- get_best_isotwas(list(a, b), G = 2)

    expect_equal(got$r2_comparison["t1", "A"], 0)
})

test_that("CHARACTERIZATION: a model selecting no SNPs is scored 0, not -Inf", {

    # Current behaviour (utils.R:69). Because calc.r2 passes negative adjusted
    # R2 through, an EMPTY model scored 0 outranks a real model scoring below
    # zero, and is then emitted as that transcript's weights. Scoring an empty
    # model -Inf would make the real model win instead.
    empty <- mk_model("A", c(0.40), n_snps = 0L)
    real <- mk_model("B", c(-0.01), n_snps = 30L)

    got <- get_best_isotwas(list(empty, real), G = 1)

    expect_equal(got$r2_comparison["t1", "A"], 0)
    expect_identical(got$best_methods$best_method, "A")
    expect_identical(nrow(got$transcripts[["t1"]]$weights), 0L)
})

test_that("CHARACTERIZATION: a nested ensemble's leaf label is overwritten", {

    # Current behaviour (utils.R:88). compute_isotwas builds its `univariate`
    # entry by calling get_best() a second time, so that entry's own $method is
    # the literal "best". At the outer level every transcript won by the
    # ensemble is relabelled "best", destroying which of enet/BLUP/SuSiE won.
    flat <- mk_model("spls", c(0.10, 0.90))
    ens <- mk_ensemble(c(0.50, 0.40),
                       leaves = c("univariate_blup", "univariate_susie"))

    got <- get_best_isotwas(list(flat, ens), G = 2)

    expect_identical(got$best_methods$best_method, c("best", "spls"))
    expect_identical(got$transcripts[["t1"]]$method, "best")
})
