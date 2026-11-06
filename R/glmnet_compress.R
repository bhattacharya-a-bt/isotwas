#' Exact glmnet fits from compressed rows
#'
#' A gaussian or mgaussian elastic-net path depends on the data only through
#' the centred cross-products X'X, X'Y and Y'Y and the number of rows n: glmnet
#' centres x and y (intercept = TRUE), optionally divides each column of x by
#' its population standard deviation, and minimises
#' (1/2n) RSS + lambda * penalty. Any k rows with the same cross-products,
#' scaled by sqrt(k / n) so that (1/2k) RSS equals the original (1/2n) RSS,
#' give glmnet the same problem. When n is well above p, k = rank(X) + 1 rows
#' per response (rank(X) + q for mgaussian) make every pass of the coordinate
#' descent cheaper, and the path is the same to rounding error.
#'
#' The rows come from a pivoted Cholesky factor of the centred X'X: with
#' C'C = X'X (r = rank(X) rows, columns unpivoted) and Q'Y = C^-T X'Y on the r
#' pivot columns, response j gets the rows [C, (Q'Y)[, j]] plus one row
#' [0, e_j], where e_j^2 = y_j'y_j - |(Q'Y)[, j]|^2 is the residual sum of
#' squares of y_j on X; mgaussian gets [C, Q'Y] over [0, R_e], with R_e'R_e
#' the residual cross-product matrix. Standardising divides the columns of C
#' by the standard deviations, which leaves Q'Y unchanged.
#'
#' The cross-products of any subset of rows (a CV training fold) follow from
#' those of the full data minus those of the dropped rows, so each fold costs
#' a product over the held-out rows plus a p^3/3 factorisation, shared by
#' every response, alpha and method that uses that fold.
#'
#' glmnet is then called with intercept = FALSE, standardize = FALSE and the
#' lambda.min.ratio its own default would pick for the original n and p (the
#' default depends on the number of rows). The fit is repackaged on the
#' original scale (beta / sd, intercept ym - xm' beta, nobs, nulldev), so
#' coef() and predict() work on it as on the fit to the original rows.
#'
#' Compression is on with options(isotwas.glmnet_compress = TRUE), and by
#' default inside a compute_isotwas() run with shared folds, where each
#' design serves several methods (see .glmnet_memo_start()); it is off
#' otherwise, and everywhere with options(isotwas.glmnet_compress = FALSE).
#' When on, it is used for every design that removes enough rows: p + q + 1
#' at most getOption('isotwas.glmnet_compress_ratio', 0.75) times the rows
#' of the subset (.glmnet_compress_on()). Whether that is faster depends on
#' the machine: building a design is BLAS work (cross-products and a p^3/3
#' factorisation), cheap with an optimised BLAS and some 30 times dearer
#' with R's reference BLAS, where a design that serves only one or two fits
#' can cost more than it saves. The choice depends only on the option, the
#' call and the sizes, never on timings, so a rerun takes the same route
#' and gives identical results; the two routes agree to rounding error.
#'
#' @noRd
NULL

#' Use compressed rows for a design?
#'
#' options(isotwas.glmnet_compress) TRUE or FALSE decides; unset, the
#' current compute_isotwas run does (.glmnet_memo_start(compress)), and
#' outside one compression is off.
#'
#' @param n rows of the subset the design is for
#' @param p SNPs
#' @param q responses per fit (rows added per fit)
#' @noRd
.glmnet_compress_on <- function(n, p, q = 1L) {
    mode <- getOption("isotwas.glmnet_compress")
    on <- if (is.null(mode)) isTRUE(.isotwas_cache$memo$compress)
          else isTRUE(mode)
    on && p + q + 1 <= getOption("isotwas.glmnet_compress_ratio", 0.75) * n
}

#' Cross-products of the full data, for downdating to subsets
#'
#' X and Y are shifted by their full-data column means first, so the
#' subset means are small and centring a subset by its own means loses
#' little precision.
#'
#' @return list with the shifted X0 and Y0, their raw cross-products and
#'   column sums, and the shifts
#' @noRd
.gram_full <- function(X, Y) {
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    mx <- colMeans(X)
    my <- colMeans(Y)
    X0 <- sweep(X, 2, mx)
    Y0 <- sweep(Y, 2, my)
    list(X0 = X0, Y0 = Y0, mx = mx, my = my, n = nrow(X),
         sxx = crossprod(X0), sxy = crossprod(X0, Y0), syy = crossprod(Y0),
         sx = colSums(X0), sy = colSums(Y0))
}

#' Centred cross-products of the rows that remain after dropping some
#'
#' @param g output of .gram_full()
#' @param drop integer or logical index of the rows to leave out
#' @return list(n, xm, ym, xx, xy, yy), cross-products centred at the
#'   subset's own means
#' @noRd
.gram_subset <- function(g, drop = NULL) {
    sxx <- g$sxx; sxy <- g$sxy; syy <- g$syy; sx <- g$sx; sy <- g$sy
    n <- g$n
    if (length(drop) > 0) {
        if (is.logical(drop)) drop <- which(drop)
        Xd <- g$X0[drop, , drop = FALSE]
        Yd <- g$Y0[drop, , drop = FALSE]
        sxx <- sxx - crossprod(Xd)
        sxy <- sxy - crossprod(Xd, Yd)
        syy <- syy - crossprod(Yd)
        sx <- sx - colSums(Xd)
        sy <- sy - colSums(Yd)
        n <- n - length(drop)
    }
    list(n = n, xm = g$mx + sx / n, ym = g$my + sy / n,
         xx = sxx - tcrossprod(sx) / n,
         xy = sxy - tcrossprod(sx, sy) / n,
         yy = syy - tcrossprod(sy) / n)
}

#' Factor centred cross-products into compressed rows
#'
#' @param s output of .gram_subset()
#' @return list(n, p, xm, ym, xs, R, QtY, yy, enorm, Re); R is rank(X) x p
#'   on the centred, unscaled scale
#' @noRd
.glmnet_design <- function(s) {
    p <- ncol(s$xx)
    C <- suppressWarnings(chol(s$xx, pivot = TRUE))
    r <- attr(C, "rank")
    piv <- attr(C, "pivot")
    C <- C[seq_len(r), , drop = FALSE]
    QtY <- backsolve(C[, seq_len(r), drop = FALSE],
                     s$xy[piv[seq_len(r)], , drop = FALSE],
                     transpose = TRUE)
    E <- s$yy - crossprod(QtY)
    Ce <- suppressWarnings(chol(E, pivot = TRUE))
    re <- attr(Ce, "rank")
    if (re < ncol(E)) Ce[-seq_len(re), ] <- 0
    xs <- sqrt(diag(s$xx) / s$n)
    xs[xs == 0] <- 1
    R <- C[, order(piv), drop = FALSE]
    colnames(R) <- colnames(s$xx)
    list(n = s$n, p = p, xm = s$xm, ym = s$ym, xs = xs, R = R,
         QtY = QtY, yy = s$yy, enorm = sqrt(pmax(diag(E), 0)),
         Re = Ce[, order(attr(Ce, "pivot")), drop = FALSE])
}

#' glmnet on compressed rows, repackaged on the original scale
#'
#' @param d output of .glmnet_design()
#' @param j response column for family 'gaussian'; ignored for 'mgaussian'
#' @param lambda NULL for glmnet's own path, or user lambdas as in glmnet()
#' @param ... further glmnet() arguments (e.g. lower.limits)
#' @return an elnet or mrelnet object equal, to rounding, to
#'   glmnet(X_subset, Y_subset, family, alpha = alpha,
#'   standardize = standardize, lambda = lambda, ...)
#' @noRd
.glmnet_compressed <- function(d, j = 1L, family = c("gaussian", "mgaussian"),
                               alpha = 1, standardize = TRUE, lambda = NULL,
                               ...) {
    family <- match.arg(family)
    sc <- if (standardize) d$xs else rep(1, d$p)
    Rs <- sweep(d$R, 2, sc, "/")
    r <- nrow(Rs)
    if (family == "gaussian") {
        k <- r + 1
        xr <- rbind(Rs, 0)
        yr <- c(d$QtY[, j], d$enorm[j])
    } else {
        q <- ncol(d$QtY)
        k <- r + q
        xr <- rbind(Rs, matrix(0, q, d$p))
        yr <- rbind(d$QtY, d$Re)
        colnames(yr) <- colnames(d$yy)
    }
    w <- sqrt(k / d$n)
    xr <- xr * w
    yr <- yr * w
    colnames(xr) <- colnames(d$R)
    lmr <- if (d$n < d$p) 1e-2 else 1e-4
    fit <- glmnet::glmnet(xr, yr, family = family, alpha = alpha,
                          standardize = FALSE, intercept = FALSE,
                          lambda.min.ratio = lmr, lambda = lambda, ...)
    if (family == "gaussian") {
        fit$beta <- Matrix::Diagonal(x = 1 / sc) %*% fit$beta
        dimnames(fit$beta) <- list(colnames(d$R), names(fit$a0))
        a0 <- d$ym[j] - as.numeric(crossprod(d$xm, as.matrix(fit$beta)))
        names(a0) <- names(fit$a0)
        fit$a0 <- a0
        fit$nulldev <- d$yy[j, j]
    } else {
        dn <- dimnames(fit$beta[[1]])
        for (m in seq_along(fit$beta)) {
            b <- Matrix::Diagonal(x = 1 / sc) %*% fit$beta[[m]]
            dimnames(b) <- dn
            fit$beta[[m]] <- b
            fit$a0[m, ] <- d$ym[m] - as.numeric(crossprod(d$xm, as.matrix(b)))
        }
        fit$nulldev <- sum(diag(d$yy))
    }
    fit$nobs <- d$n
    fit
}

.isotwas_cache <- new.env(parent = emptyenv())

#' Designs shared across methods within one compute_isotwas run
#'
#' compute_isotwas fits several glmnet-based methods to the same X and Y on
#' the same (shared) folds. While .glmnet_memo_start() is in effect,
#' .gram_full_memo() returns the cross-products already formed for an
#' identical X and Y, and .gram_design() the design already built for the
#' same rows, so each is built once per run. The memo holds one X and Y at
#' a time and is cleared by .glmnet_memo_end().
#'
#' With compress = TRUE (compute_isotwas passes it when its folds are
#' shared), compression is on for the run unless
#' options(isotwas.glmnet_compress) says otherwise (.glmnet_compress_on()).
#' Each design then serves every method, alpha and transcript that uses its
#' rows, which makes it worth building even with R's reference BLAS.
#' @param compress logical, compress by default during the run
#' @noRd
.glmnet_memo_start <- function(compress = FALSE) {
    .isotwas_cache$memo <- new.env(parent = emptyenv())
    .isotwas_cache$memo$compress <- compress
    invisible(NULL)
}

#' @rdname dot-glmnet_memo_start
#' @noRd
.glmnet_memo_end <- function() {
    .isotwas_cache$memo <- NULL
    invisible(NULL)
}

.memo_key <- function(drop) {
    if (is.logical(drop)) drop <- which(drop)
    paste0("d", paste(sort(drop), collapse = ","))
}

.memo_matches <- function(m, X, Y) {
    !is.null(m) && !is.null(m$X) &&
        identical(m$X, X) && identical(m$Y, unname(as.matrix(Y)))
}

#' .gram_full(), from the memo when it holds the same X and Y
#' @noRd
.gram_full_memo <- function(X, Y) {
    m <- .isotwas_cache$memo
    if (.memo_matches(m, X, Y)) return(m$gram)
    g <- .gram_full(X, Y)
    if (!is.null(m)) {
        m$X <- X
        m$Y <- unname(as.matrix(Y))
        m$gram <- g
        m$designs <- list()
    }
    g
}

#' .glmnet_design(.gram_subset(g, drop)), from the memo when built already
#' @param keep logical, store a newly built design in the memo; FALSE for
#'   rows no other fit will use (folds drawn for one transcript)
#' @noRd
.gram_design <- function(g, drop = integer(0), keep = TRUE) {
    m <- .isotwas_cache$memo
    key <- .memo_key(drop)
    shared <- !is.null(m) && identical(m$gram, g)
    if (shared && !is.null(m$designs[[key]])) return(m$designs[[key]])
    d <- .glmnet_design(.gram_subset(g, drop))
    if (shared && keep) m$designs[[key]] <- d
    d
}


#' Compressed designs for the full data and each CV training fold
#'
#' @param X,Y full data (Y n x q)
#' @param foldid fold labels 1..K
#' @param g,full optional .gram_full(X, Y) and full-data design to reuse
#'   across several fold sets
#' @param keep logical, store the fold designs in the run's memo
#' @return list(full, folds): .glmnet_design() outputs
#' @noRd
.glmnet_cv_designs <- function(X, Y, foldid, g = NULL, full = NULL,
                               keep = TRUE) {
    if (is.null(g)) g <- .gram_full_memo(X, Y)
    if (is.null(full)) full <- .gram_design(g)
    list(full = full,
         folds = lapply(seq_len(max(foldid)), function(i)
             .gram_design(g, which(foldid == i), keep = keep)))
}

#' cv.glmnet for gaussian and mgaussian fits, on compressed rows
#'
#' Reproduces glmnet::cv.glmnet(X, y, foldid = foldid, family = family,
#' alpha = alpha, standardize = standardize, keep = TRUE) with the default
#' type.measure (mean squared error), grouped = TRUE and alignment = 'lambda'
#' (glmnet 4.1 to 5.0): a fit to all rows, a fit with its own lambda path to
#' each training fold, held-out predictions at the full-data lambdas
#' (interpolated, as predict.glmnet does), the fold-size-weighted mean of the
#' per-fold mean squared errors, and lambda.min as the largest lambda at the
#' minimum. The caller draws foldid, as cv.glmnet's
#' sample(rep(seq(nfolds), length = N)) would.
#'
#' Fits use compressed rows when compression is on (see above) and the
#' original rows
#' otherwise; the returned object has the class and fields of cv.glmnet's,
#' so coef(), predict() and fit.preval work as before.
#'
#' @param X n x p matrix
#' @param Y response: vector (gaussian) or n x q matrix (mgaussian)
#' @param foldid fold labels 1..K
#' @param designs optional list(full, folds) of .glmnet_design() outputs to
#'   reuse, built from these X and Y columns and foldid; a NULL entry means
#'   that fit uses the original rows. Without it, .cv_glmnet builds the
#'   designs when compression is on.
#' @param j column of designs' responses that Y is, for family 'gaussian'
#' @param n.cores folds fitted in parallel with mclapply when > 1
#' @noRd
.cv_glmnet <- function(X, Y, foldid, family = c("gaussian", "mgaussian"),
                       alpha = 1, standardize = TRUE, designs = NULL, j = 1L,
                       n.cores = 1L) {
    family <- match.arg(family)
    X <- as.matrix(X)
    y <- drop(Y)
    N <- nrow(X)
    nfolds <- max(foldid)
    if (nfolds < 3)
        stop("nfolds must be bigger than 3; nfolds=10 recommended")
    q <- if (family == "mgaussian") NCOL(y) else 1L
    # one design for all rows and one per fold
    if (is.null(designs) &&
        .glmnet_compress_on(N - max(tabulate(foldid)), ncol(X), q)) {
        g <- .gram_full(X, as.matrix(Y))
        designs <- list(full = .glmnet_design(.gram_subset(g)),
                        folds = lapply(seq_len(nfolds), function(i)
                            .glmnet_design(.gram_subset(g, foldid == i))))
        j <- 1L
    }

    fit_one <- function(i) {
        d <- if (is.null(designs)) NULL
             else if (i == 0) designs$full else designs$folds[[i]]
        if (!is.null(d)) {
            .glmnet_compressed(d, j, family, alpha, standardize)
        } else if (i == 0) {
            glmnet::glmnet(X, y, family = family, alpha = alpha,
                           standardize = standardize)
        } else {
            keep <- foldid != i
            y_sub <- if (length(dim(y)) > 1) y[keep, ] else y[keep]
            glmnet::glmnet(X[keep, , drop = FALSE], y_sub, family = family,
                           alpha = alpha, standardize = standardize)
        }
    }

    full <- fit_one(0)
    lambda <- full$lambda
    nz <- sapply(stats::predict(full, type = "nonzero"), length)
    outlist <- if (n.cores > 1) {
        parallel::mclapply(seq_len(nfolds), fit_one, mc.cores = n.cores)
    } else {
        lapply(seq_len(nfolds), fit_one)
    }

    # buildPredmat: held-out predictions at the full-data lambdas
    nlambda <- length(lambda)
    sn <- paste("s", seq(0, length = nlambda), sep = "")
    if (family == "gaussian") {
        predmat <- matrix(NA, N, nlambda)
        for (i in seq_len(nfolds)) {
            which <- foldid == i
            preds <- stats::predict(outlist[[i]], X[which, , drop = FALSE],
                                  s = lambda)
            nlami <- min(ncol(preds), nlambda)
            predmat[which, seq(nlami)] <- preds[, seq(nlami)]
            if (nlami < nlambda)
                predmat[which, seq(from = nlami, to = nlambda)] <- preds[, nlami]
        }
        dimnames(predmat) <- list(rownames(X), sn)
        cvraw <- (y - predmat)^2
    } else {
        nc <- ncol(y)
        predmat <- array(NA, c(N, nc, nlambda))
        for (i in seq_len(nfolds)) {
            which <- foldid == i
            preds <- stats::predict(outlist[[i]], X[which, , drop = FALSE],
                                  s = lambda)
            nlami <- min(dim(preds)[3], nlambda)
            predmat[which, , seq(nlami)] <- preds[, , seq(nlami)]
            if (nlami < nlambda)
                predmat[which, , seq(from = nlami, to = nlambda)] <- preds[, , nlami]
        }
        dimnames(predmat) <- list(rownames(X), dimnames(preds)[[2]], sn)
        cvraw <- apply((array(y, dim(predmat)) - predmat)^2, c(1, 3), sum)
    }

    # cvcompute and cvstats (grouped = TRUE, unit weights)
    weights <- rep(1, N)
    grouped <- !(N / nfolds < 3)
    if (!grouped)
        warning("Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold",
                call. = FALSE)
    if (grouped) {
        wisum <- tapply(weights, foldid, sum)
        outmat <- matrix(NA, nfolds, ncol(cvraw))
        cvraw[is.infinite(cvraw)] <- NA
        for (i in seq_len(nfolds)) {
            outmat[i, ] <- apply(cvraw[foldid == i, , drop = FALSE], 2,
                                 stats::weighted.mean, w = weights[foldid == i],
                                 na.rm = TRUE)
        }
        cvraw <- outmat
        weights <- wisum
        Nobs <- rep(nfolds, ncol(cvraw))
    } else {
        Nobs <- N - apply(is.na(cvraw), 2, sum)
    }
    cvm <- apply(cvraw, 2, stats::weighted.mean, w = weights, na.rm = TRUE)
    cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, stats::weighted.mean,
                       w = weights, na.rm = TRUE) / (Nobs - 1))
    nas <- is.na(cvsd)
    if (any(nas)) {
        lambda <- lambda[!nas]
        cvm <- cvm[!nas]
        cvsd <- cvsd[!nas]
        nz <- nz[!nas]
    }
    cvname <- "Mean-Squared Error"
    names(cvname) <- "mse"

    # getOptcv.glmnet
    lambda.min <- max(lambda[cvm <= min(cvm, na.rm = TRUE)], na.rm = TRUE)
    idmin <- match(lambda.min, lambda)
    lambda.1se <- max(lambda[cvm <= (cvm + cvsd)[idmin]], na.rm = TRUE)
    id1se <- match(lambda.1se, lambda)

    obj <- list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvup = cvm + cvsd,
                cvlo = cvm - cvsd, nzero = nz, call = match.call(),
                name = cvname, glmnet.fit = full, fit.preval = predmat,
                foldid = foldid, lambda.min = lambda.min,
                lambda.1se = lambda.1se,
                index = matrix(c(idmin, id1se), 2, 1,
                               dimnames = list(c("min", "1se"), "Lambda")))
    class(obj) <- "cv.glmnet"
    obj
}
