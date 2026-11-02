# A model must report the name it appears under in the comparison table.
#
# These disagreed: multivariate_elasticnet() labelled itself
# "multivariate_enet" against a column named `multi_enet`, multivariate_sgl()
# "sparse_group_lasso" against `sgl`, multivariate_mtlasso() "mtlasso_l21"
# against `mtlasso`, and multivariate_graph_reg() "graph_regularized" against
# `graph_reg`. Every consumer of a fitted panel therefore needed its own alias
# table to work out which method produced a given transcript's weights.

test_that("the method list is a set, and matches the documented argument", {

    m <- isotwas_methods()

    expect_type(m, "character")
    expect_false(anyDuplicated(m) > 0)
    expect_true(all(nzchar(m)))

    # The `method` argument of compute_isotwas() is documented as taking these.
    expect_setequal(m, c('mrce_lasso', 'multi_enet', 'joinet', 'spls', 'sgl',
                         'mtlasso', 'stacking', 'graph_reg', 'univariate'))
})

test_that("multivariate_elasticnet labels itself 'multi_enet'", {

    skip_if_not_installed("glmnet")
    skip_if_not_installed("tibble")

    d <- withr::with_seed(51, mk_xy(n = 60, p = 15, q = 3))

    mod <- multivariate_elasticnet(X = d$X, Y = d$Y, Omega = NULL,
                                   nfolds = 3, verbose = FALSE,
                                   par = FALSE, seed = 1)

    expect_identical(mod$method, "multi_enet")
    expect_true(mod$method %in% isotwas_methods())
})

test_that("multivariate_sgl labels itself 'sgl'", {

    skip_if_not_installed("tibble")

    d <- withr::with_seed(52, mk_xy(n = 50, p = 12, q = 3))

    mod <- multivariate_sgl(X = d$X, Y = d$Y,
                            nfolds = 3, nlambda = 5, verbose = FALSE, seed = 1)

    expect_identical(mod$method, "sgl")
    expect_true(mod$method %in% isotwas_methods())
})

test_that("every fitting function's label is one of the known methods", {

    # Guards against a new method being added with a name that does not appear
    # as a comparison column. Reads the sources rather than fitting, so it
    # covers the methods whose dependencies are not installed here.
    r_dir <- system.file("R", package = "isotwas")
    if (!nzchar(r_dir) || !dir.exists(r_dir)) r_dir <- "../../R"
    skip_if_not(dir.exists(r_dir), "package sources not available")

    files <- list.files(r_dir, pattern = "[.]R$", full.names = TRUE)
    lines <- unlist(lapply(files, readLines, warn = FALSE))

    # `method = "..."` as passed to create_isotwas_model().
    hits <- regmatches(lines, regexpr('method = "[^"]+"', lines))
    labels <- unique(gsub('method = "([^"]+)"', "\\1", hits))

    # The univariate leaves are members of the `univariate` ensemble rather
    # than comparison columns in their own right, and "best" is the assembled
    # model. Everything else must be a column.
    leaves <- c("univariate_enet", "univariate_blup", "univariate_susie",
                "best", "mean", "replicates")
    labels <- setdiff(labels, leaves)

    expect_true(all(labels %in% isotwas_methods()),
                info = paste("unknown:",
                             paste(setdiff(labels, isotwas_methods()),
                                   collapse = ", ")))

    # A computed label sidesteps the check above. The RMTL backend of
    # multivariate_mtlasso() built its name as
    # paste0("mtlasso_", tolower(regularization)), producing mtlasso_l21 and
    # friends against a column named `mtlasso`; the variant is recorded in
    # $regularization instead. Labels must be literals.
    computed <- grep("method = paste|method = sprintf|method = if",
                     lines, value = TRUE)
    expect_identical(length(computed), 0L,
                     info = paste(computed, collapse = " | "))
})

test_that("mtlasso records its regularization on every transcript", {

    skip_if_not_installed("tibble")

    d <- withr::with_seed(53, mk_xy(n = 50, p = 12, q = 3))

    # Selection copies transcript models out of their parent, so the variant
    # must live on each transcript to survive into the assembled model.
    mod <- multivariate_mtlasso(X = d$X, Y = d$Y, nfolds = 3, nlambda = 5,
                                seed = 1)
    expect_identical(mod$method, "mtlasso")
    expect_true(all(vapply(mod$transcripts, `[[`, "", "regularization") ==
                    "L21"))

    best <- get_best_isotwas(list(mod), G = 3)
    expect_true(all(vapply(best$transcripts, `[[`, "", "regularization") ==
                    "L21"))

    skip_if_not_installed("RMTL")
    rmtl <- multivariate_mtlasso(X = d$X, Y = d$Y, regularization = "Trace",
                                 nfolds = 3, seed = 1)
    expect_identical(rmtl$method, "mtlasso")
    expect_true(all(vapply(rmtl$transcripts, `[[`, "", "regularization") ==
                    "Trace"))
})
