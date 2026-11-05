#' spls::cv.spls and spls::spls as multivariate_spls calls them
#'
#' multivariate_spls runs spls::cv.spls(x, y, fold, K = 1:max_K, eta =
#' seq(.1, .9, .1), scale.x = FALSE, scale.y = FALSE) and then spls::spls()
#' with the selected K and eta, once on all samples and once per fold. These
#' functions reproduce spls 2.3-2 (select = "pls2", fit = "simpls", kappa =
#' 0.5) with the same arithmetic in the same order, so their results are the
#' same to the last bit. They are faster because they
#' \itemize{
#'   \item call pls::simpls.fit directly, the routine pls::plsr dispatches to,
#'     instead of going through a model frame on every fit;
#'   \item centre each training set, and form its X'Y and the spls.dv matrix
#'     Z Z' of the first component, once for all eta rather than once per eta;
#'   \item work on numeric matrices rather than data frames.
#' }
#'
#' @noRd
NULL

#' spls:::ust
#' @noRd
.spls_ust <- function(b, eta) {
    b.ust <- matrix(0, length(b), 1)
    if (eta < 1) {
        valb <- abs(b) - eta * max(abs(b))
        b.ust[valb >= 0] <- valb[valb >= 0] * (sign(b))[valb >= 0]
    }
    b.ust
}

#' spls:::spls.dv with kappa = 0.5; Zn and M may be passed in precomputed
#' @noRd
.spls_dv <- function(Z, eta, eps, maxstep, Zn = NULL, M = NULL) {
    p <- nrow(Z)
    q <- ncol(Z)
    if (is.null(Zn)) {
        Znorm1 <- stats::median(abs(Z))
        Zn <- Z / Znorm1
    }
    if (q == 1) {
        return(.spls_ust(Zn, eta))
    }
    # M = Zn Zn' is p x p; M v is formed as Zn (Zn' v), O(pq) rather than
    # O(p^2), which changes only the rounding
    Mv <- function(v) Zn %*% crossprod(Zn, v)
    dis <- 10
    i <- 1
    c <- matrix(10, p, 1)
    c.old <- c
    while (dis > eps & i <= maxstep) {
        mcsvd <- svd(Mv(c))
        a <- mcsvd$u %*% t(mcsvd$v)
        c <- .spls_ust(Mv(a), eta)
        dis <- max(abs(c - c.old))
        c.old <- c
        i <- i + 1
    }
    c
}

#' Centre a training set as spls::spls does, and precompute the first
#' component's direction-vector inputs, which do not depend on eta
#' @noRd
.spls_prep <- function(x, y) {
    x <- as.matrix(x)
    y <- as.matrix(y)
    n <- nrow(x)
    one <- matrix(1, 1, n)
    mu <- one %*% y / n
    y <- scale(y, drop(mu), FALSE)
    meanx <- drop(one %*% x) / n
    x <- scale(x, meanx, FALSE)
    Z <- t(x) %*% y
    Zn <- Z / stats::median(abs(Z))
    M <- NULL
    list(x = x, y = y, mu = mu, meanx = meanx, Z = Z, Zn = Zn, M = M)
}

#' spls::spls(x, y, K, eta, scale.x = FALSE, scale.y = FALSE) on a prepared
#' training set
#'
#' @return list(betahat, betamat, A, mu, meanx), as the fields of the spls
#'   object of the same names
#' @noRd
.spls_fit <- function(prep, K, eta, eps = 1e-04, maxstep = 100) {
    x <- prep$x
    y <- prep$y
    n <- nrow(x)
    p <- ncol(x)
    q <- ncol(y)
    ip <- c(1:p)
    if (eta < 0 || eta >= 1) stop("eta should be between 0 and 1!")
    if (K > p) stop("K cannot exceed the number of predictors! Pick up smaller K!")
    if (K >= n) stop("K cannot exceed the sample size! Pick up smaller K!")

    betahat <- matrix(0, p, q)
    betamat <- list()
    y1 <- y
    for (k in 1:K) {
        what <- if (k == 1) {
            .spls_dv(prep$Z, eta, eps, maxstep, Zn = prep$Zn, M = prep$M)
        } else {
            .spls_dv(t(x) %*% y1, eta, eps, maxstep)
        }
        A <- unique(ip[what != 0 | betahat[, 1] != 0])
        xA <- x[, A, drop = FALSE]
        ncomp <- min(k, length(A))
        if (ncomp > min(n - 1, length(A))) stop("Invalid number of components, ncomp")
        fit <- pls::simpls.fit(xA, y, ncomp, center = TRUE, stripped = TRUE)
        betahat <- matrix(0, p, q)
        betahat[A, ] <- matrix(fit$coefficients[, , ncomp], length(A), q)
        betamat[[k]] <- betahat
        y1 <- y - x %*% betahat
    }
    list(betahat = betahat, betamat = betamat, A = A,
         mu = prep$mu, meanx = prep$meanx)
}

#' spls::cv.spls(x, y, fold, K, eta, scale.x = FALSE, scale.y = FALSE,
#' plot.it = FALSE), without its progress printing
#'
#' Draws its folds with the same single call to sample() as cv.spls.
#' @return list(mspemat, eta.opt, K.opt)
#' @noRd
.cv_spls <- function(x, y, fold, K, eta) {
    x <- as.matrix(x)
    y <- as.matrix(y)
    n <- nrow(x)
    if (max(K) > ncol(x)) stop("K cannot exceed the number of predictors! Pick up smaller K!")
    if (max(K) >= n) stop("K cannot exceed the sample size! Pick up smaller K!")
    foldi <- split(sample(1:n), rep(1:fold, length = n))
    mspemat <- matrix(0, length(eta), length(K))
    preps <- lapply(seq_len(fold), function(j) {
        omit <- foldi[[j]]
        .spls_prep(x[-omit, , drop = FALSE], y[-omit, , drop = FALSE])
    })
    for (i in 1:length(eta)) {
        mspemati <- matrix(0, fold, length(K))
        for (j in 1:fold) {
            omit <- foldi[[j]]
            object <- .spls_fit(preps[[j]], K = max(K), eta = eta[i])
            newx <- x[omit, , drop = FALSE]
            newx <- scale(newx, object$meanx, rep(1, ncol(x)))
            betamat <- object$betamat
            for (k in K) {
                pred <- newx %*% betamat[[k]] + matrix(1, nrow(newx), 1) %*% object$mu
                mspemati[j, (k - min(K) + 1)] <- mean(apply((y[omit, ] - pred)^2, 2, mean))
            }
        }
        mspemat[i, ] <- apply(mspemati, 2, mean)
    }
    minpmse <- min(mspemat)
    mspecol <- apply(mspemat, 2, min)
    msperow <- apply(mspemat, 1, min)
    K.opt <- min(K[mspecol == minpmse])
    eta.opt <- max(eta[msperow == minpmse])
    rownames(mspemat) <- paste("eta=", eta)
    colnames(mspemat) <- paste("K =", K)
    list(mspemat = mspemat, eta.opt = eta.opt, K.opt = K.opt)
}
