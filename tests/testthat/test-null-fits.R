# lm(obs ~ pred) refits the slope, so an anti-correlated cross-validated
# prediction used to score a positive R2 and could win selection.

test_that("an anti-correlated prediction scores 0", {

    y <- withr::with_seed(71, stats::rnorm(60))
    good <- y + withr::with_seed(72, stats::rnorm(60))
    flipped <- -good

    expect_gt(.score_r2(y, good)$R2, 0)
    expect_identical(.score_r2(y, flipped)$R2, 0)
    expect_identical(.score_r2(y, rep(1, 60))$R2, 0)

    Y <- cbind(y, y); P <- cbind(good, flipped)
    expect_gt(calc.r2(1, Y, P)$R2, 0)
    expect_identical(calc.r2(2, Y, P)$R2, 0)
})

test_that(".drop_null_fits empties null fits and keeps real ones", {

    y <- withr::with_seed(73, stats::rnorm(60))
    pred <- y + withr::with_seed(74, stats::rnorm(60))
    Y <- cbind(t1 = y, t2 = y, t3 = y)
    w <- function(v) data.frame(SNP = paste0("rs", seq_along(v)), Weight = v)

    mod <- create_isotwas_model(
        method = "multi_enet",
        transcripts = list(
            t1 = create_transcript_model("t1", w(c(0.2, -0.1)), 0.4, 0.01, pred),
            t2 = create_transcript_model("t2", w(c(0.2, -0.1)), 0.3, 0.01, -pred),
            t3 = create_transcript_model("t3", w(c(1e-38, 2e-38)), 0.3, 0.01, pred)),
        n_samples = 60, n_snps = 2)

    out <- .drop_null_fits(mod, Y)

    expect_identical(nrow(out$transcripts$t1$weights), 2L)
    expect_identical(out$transcripts$t1$r2, 0.4)
    for (tx in c("t2", "t3")) {
        expect_identical(nrow(out$transcripts[[tx]]$weights), 0L)
        expect_identical(out$transcripts[[tx]]$r2, 0)
    }
    expect_identical(out$summary$n_snps_selected, c(2L, 0L, 0L))
    expect_identical(out$summary$r2, c(0.4, 0, 0))
})
