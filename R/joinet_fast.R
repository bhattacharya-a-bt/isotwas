#' joinet 1.0.0 for gaussian outputs, as multivariate_joinet calls it
#'
#' multivariate_joinet ran joinet::cv.joinet(Y, X, nfolds.ext, nfolds.int,
#' alpha.base = alpha.meta = alpha, cvpred = TRUE) for the cross-validated
#' predictions and then joinet::joinet(Y, X, nfolds = 5) for the weights.
#' The functions here reproduce both for family = 'gaussian' with no
#' missing values:
#' \itemize{
#'   \item .joinet_fit() is joinet::joinet(): the sign matrix from
#'     Spearman tests (p <= 0.05), one glmnet path per output on all rows
#'     (alpha.base, standardised), the same fit on each training fold
#'     predicted on its held-out rows at the full-data lambdas, lambda.min
#'     from the pooled mean squared error, and a meta cv.glmnet per output
#'     on the held-out predictions, with sign limits, alpha.meta and the
#'     same folds.
#'   \item .cv_joinet() is cv.joinet(cvpred = TRUE)$cvpred: external folds,
#'     constant SNPs dropped within each, internal folds per external fold,
#'     and .joinet_fit() on each external training set. cv.joinet also
#'     fitted an unused mgaussian glmnet in every external fold; that is
#'     skipped (it draws no random numbers).
#'   \item .joinet_folds() is palasso:::.folds(), the only random draw.
#' }
#' The base paths use compressed rows when compression is on
#' (glmnet_compress.R),
#' with one set of cross-products per call shared by all outputs; the meta
#' models (q columns) are fitted by cv.glmnet as before.
#'
#' @noRd
NULL

#' palasso:::.folds (palasso 1.0.0) without the foldid and Surv branches
#' @noRd
.joinet_folds <- function(y, nfolds) {
    if (all(y %in% c(0, 1))) {
        foldid <- rep(x = NA, times = length(y))
        foldid[y == 0] <- sample(x = rep(x = seq_len(nfolds),
                                         length.out = sum(y == 0)))
        foldid[y == 1] <- sample(x = rep(x = seq_len(nfolds),
                                         length.out = sum(y == 1)))
    } else {
        foldid <- sample(x = rep(x = seq_len(nfolds), length.out = length(y)))
    }
    foldid
}

#' joinet::joinet for gaussian outputs with given folds
#'
#' @param X n x p matrix
#' @param Y n x q matrix without missing values
#' @param foldid fold labels 1..K
#' @return list(base, meta): base[[i]] holds glmnet.fit, lambda, cvm and
#'   lambda.min; meta[[i]] is the meta cv.glmnet
#' @noRd
.joinet_fit <- function(X, Y, foldid, alpha.base = 1, alpha.meta = 1) {
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    n <- nrow(Y)
    q <- ncol(Y)
    nfolds <- length(unique(foldid))

    sign <- matrix(0, nrow = q, ncol = q)
    diag(sign) <- 1
    for (i in seq_len(q)) {
        for (j in seq(from = i, to = q, by = 1)) {
            cor <- stats::cor.test(Y[, i], Y[, j],
                                   use = "pairwise.complete.obs",
                                   method = 'spearman', exact = FALSE)
            if (is.na(cor$p.value) | cor$p.value > 0.05) next
            sign[i, j] <- sign[j, i] <- sign(cor$estimate)
        }
    }

    # each design (full data or one fold) serves all q outputs
    compress <- .glmnet_compress_on(n - max(tabulate(foldid)), ncol(X))
    g <- full <- dk <- NULL
    if (compress) {
        g <- .gram_full(X, Y)
        full <- .glmnet_design(.gram_subset(g))
    }
    # compressed fits take a design (full data or one training fold);
    # otherwise glmnet gets the rows as joinet passed them
    fit_base <- function(i, k = 0, d = full) {
        if (compress) {
            .glmnet_compressed(d, i, "gaussian", alpha = alpha.base)
        } else if (k == 0) {
            glmnet::glmnet(y = Y[, i], x = X, family = "gaussian",
                           alpha = alpha.base, penalty.factor = rep(1, ncol(X)))
        } else {
            keep <- foldid != k
            glmnet::glmnet(y = Y[keep, i], x = X[keep, , drop = FALSE],
                           family = "gaussian", alpha = alpha.base,
                           penalty.factor = rep(1, ncol(X)))
        }
    }

    base <- lapply(seq_len(q), function(i) {
        fit <- fit_base(i)
        list(glmnet.fit = fit, lambda = fit$lambda)
    })
    link <- lapply(base, function(b) matrix(NA, nrow = n, ncol = length(b$lambda)))
    for (k in seq_len(nfolds)) {
        which <- foldid == k
        if (compress) {
            dk <- .glmnet_design(.gram_subset(g, which))
        }
        for (i in seq_len(q)) {
            object <- fit_base(i, k, dk)
            temp <- stats::predict(object, newx = X[which, , drop = FALSE],
                                   type = "link", s = base[[i]]$lambda)
            link[[i]][which, seq_len(ncol(temp))] <- temp
        }
    }
    for (i in seq_len(q)) {
        y <- Y[, i]
        base[[i]]$cvm <- apply(link[[i]], 2, function(x) mean((x - y)^2))
        base[[i]]$lambda.min <- base[[i]]$lambda[which.min(base[[i]]$cvm)]
    }
    hat <- matrix(NA, nrow = n, ncol = q)
    for (i in seq_len(q)) {
        hat[, i] <- link[[i]][, base[[i]]$lambda == base[[i]]$lambda.min]
    }

    meta <- lapply(seq_len(q), function(i) {
        lower.limits <- rep(-Inf, times = q)
        upper.limits <- rep(Inf, times = q)
        lower.limits[sign[, i] >= 0] <- 0
        upper.limits[sign[, i] <= 0] <- 0
        glmnet::cv.glmnet(y = Y[, i], x = hat, lower.limits = lower.limits,
                          upper.limits = upper.limits, foldid = foldid,
                          family = "gaussian", type.measure = "deviance",
                          intercept = TRUE, alpha = alpha.meta)
    })
    names(base) <- names(meta) <- paste0("y", seq_len(q))
    list(base = base, meta = meta)
}

#' predict.joinet(type = 'response') for gaussian outputs
#' @noRd
.joinet_predict <- function(fit, newx) {
    newx <- as.matrix(newx)
    q <- length(fit$base)
    base <- meta <- matrix(NA, nrow = nrow(newx), ncol = q)
    for (i in seq_len(q)) {
        base[, i] <- as.numeric(stats::predict(fit$base[[i]]$glmnet.fit,
                                               newx = newx,
                                               s = fit$base[[i]]$lambda.min,
                                               type = "link"))
    }
    for (i in seq_len(q)) {
        meta[, i] <- as.numeric(stats::predict(fit$meta[[i]], newx = base,
                                               s = "lambda.min", type = "link"))
    }
    list(base = base, meta = meta)
}

#' coef.joinet: pooled intercepts and SNP weights
#' @return list(alpha, beta), beta p x q
#' @noRd
.joinet_coef <- function(fit) {
    coef <- lapply(fit$base, function(x)
        stats::coef(x$glmnet.fit, s = x$lambda.min))
    base_alpha <- sapply(coef, function(x) x[1, ])
    base_beta <- sapply(coef, function(x) x[-1, ])
    names(base_alpha) <- colnames(base_beta) <- names(fit$base)
    w <- as.matrix(do.call(cbind, lapply(fit$meta, function(x)
        stats::coef(x, s = x$lambda.min))))
    colnames(w) <- names(fit$meta)
    list(alpha = w[1, ] + base_alpha %*% w[-1, ],
         beta = base_beta %*% w[-1, ])
}

#' joinet::cv.joinet(cvpred = TRUE)$cvpred for gaussian outputs
#'
#' @param n.cores external folds fitted in parallel with mclapply when > 1
#'   (the folds are all drawn first, in cv.joinet's order)
#' @return n x q matrix of cross-validated meta-model predictions
#' @noRd
.cv_joinet <- function(X, Y, nfolds.ext = 5, nfolds.int = 10,
                       alpha.base = 1, alpha.meta = 1, n.cores = 1L) {
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    foldid.ext <- .joinet_folds(Y[, 1], nfolds.ext)
    foldid.int <- lapply(seq_len(nfolds.ext), function(i)
        .joinet_folds(Y[foldid.ext != i, 1], nfolds.int))
    run_fold <- function(i) {
        X0 <- X[foldid.ext != i, , drop = FALSE]
        X1 <- X[foldid.ext == i, , drop = FALSE]
        cond <- apply(X0, 2, function(x) stats::sd(x) != 0)
        X0 <- X0[, cond]
        X1 <- X1[, cond]
        fit <- .joinet_fit(X0, Y[foldid.ext != i, , drop = FALSE],
                           foldid.int[[i]], alpha.base, alpha.meta)
        .joinet_predict(fit, X1)$meta
    }
    folds <- if (n.cores > 1) {
        parallel::mclapply(seq_len(nfolds.ext), run_fold, mc.cores = n.cores)
    } else {
        lapply(seq_len(nfolds.ext), run_fold)
    }
    pred <- matrix(NA, nrow = nrow(Y), ncol = ncol(Y))
    for (i in seq_len(nfolds.ext)) {
        pred[foldid.ext == i, ] <- folds[[i]]
    }
    pred
}
