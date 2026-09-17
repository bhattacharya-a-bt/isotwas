# Contract tests for the isotwas_model S3 layer.
#
# These pin the shape of the objects the rest of the package builds on, so a
# refactor that renames or drops a field fails here rather than downstream in
# someone's weight table.

test_that("create_isotwas_model returns the documented shape", {

    mod <- create_isotwas_model(method = "spls",
                                transcripts = list(t1 = mk_tx(0.3),
                                                   t2 = mk_tx(0.1)),
                                n_samples = 100,
                                n_snps = 10)

    expect_s3_class(mod, "isotwas_model")
    expect_identical(mod$method, "spls")
    expect_identical(names(mod$transcripts), c("t1", "t2"))
    expect_identical(mod$n_transcripts, 2L)
    expect_identical(mod$n_samples, 100)
    expect_identical(mod$n_snps_total, 10)
})

test_that("create_isotwas_model requires a named transcript list", {

    expect_error(
        create_isotwas_model(method = "spls",
                             transcripts = list(mk_tx(0.3)),
                             n_samples = 100,
                             n_snps = 10),
        "named list")
})

test_that("the summary table carries one row per transcript", {

    mod <- create_isotwas_model(method = "spls",
                                transcripts = list(t1 = mk_tx(0.3),
                                                   t2 = mk_tx(0.1, n_snps = 4L)),
                                n_samples = 100,
                                n_snps = 10)

    expect_identical(nrow(mod$summary), 2L)
    expect_identical(names(mod$summary),
                     c("transcript", "r2", "pvalue", "n_snps_selected"))
    expect_equal(mod$summary$r2, c(0.3, 0.1))
    expect_equal(mod$summary$n_snps_selected, c(1L, 4L))
})

test_that("an empty model reports zero selected SNPs rather than failing", {

    mod <- create_isotwas_model(method = "spls",
                                transcripts = list(t1 = mk_tx(0.2, n_snps = 0L)),
                                n_samples = 100,
                                n_snps = 10)

    expect_equal(mod$summary$n_snps_selected, 0L)
})

test_that("get_weight_matrix is SNPs x transcripts and zero-fills", {

    # t1 selects rs1-rs3, t2 only rs1: the union is three rows, and t2's
    # unselected SNPs must come back as 0, not NA.
    mod <- create_isotwas_model(method = "spls",
                                transcripts = list(t1 = mk_tx(0.3, n_snps = 3L),
                                                   t2 = mk_tx(0.1, n_snps = 1L)),
                                n_samples = 100,
                                n_snps = 10)

    w <- get_weight_matrix(mod)

    expect_identical(dim(w), c(3L, 2L))
    expect_identical(colnames(w), c("t1", "t2"))
    expect_equal(w[c("rs2", "rs3"), "t2"], c(rs2 = 0, rs3 = 0))
})

test_that("get_prediction_matrix is samples x transcripts", {

    mod <- create_isotwas_model(method = "spls",
                                transcripts = list(t1 = mk_tx(0.3, n_pred = 12L),
                                                   t2 = mk_tx(0.1, n_pred = 12L)),
                                n_samples = 12,
                                n_snps = 10)

    pred <- get_prediction_matrix(mod)

    expect_identical(dim(pred), c(12L, 2L))
    expect_identical(colnames(pred), c("t1", "t2"))
})

test_that("get_weight_matrix rejects a non-model", {
    expect_error(get_weight_matrix(list(method = "spls")), "isotwas_model")
})
