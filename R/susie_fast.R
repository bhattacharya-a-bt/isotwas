#' susieR::susie as univariate_susie calls it, for many responses
#'
#' univariate_susie fits susieR::susie(X, y, estimate_residual_variance =
#' TRUE, estimate_prior_variance = FALSE, scaled_prior_variance = 0.1) with
#' the other arguments at their defaults (L = min(10, p), intercept and
#' standardize TRUE, max_iter = 100, convergence on the ELBO).
#' .susie_prep() computes what depends on X alone, once per training set;
#' .susie_fit() runs IBSS for one response.
#'
#' The algorithm is susieR 0.14.2's (the CRAN release), step for step,
#' whatever susieR version is installed, so a fit does not depend on it:
#' the same order of single-effect updates, the same log Bayes factors,
#' posterior moments, KL terms, ELBO, convergence test, and residual
#' variance update and bounds. susieR 0.16.x (GitHub) differs in four ways,
#' which this code does not follow: tol defaults to 1e-4 rather than 1e-3;
#' an ELBO decrease does not count as converged; a SNP with no variance
#' (d = 0) keeps its prior log odds rather than log posterior odds 0; and
#' null effects are trimmed at the end. What changes from susieR is how the
#' products with the standardised X are formed:
#' \itemize{
#'   \item X'R for the residual R is taken from X'y and a cached X'X b_l per
#'     effect when n >= p (a p x p product instead of two n x p ones), or
#'     from a cached X b_l per effect otherwise;
#'   \item the KL term's R'X E[b] is formed as (X'R)'E[b], and the ELBO's
#'     ||X b_l||^2 from the same caches.
#' }
#' The estimates agree with susieR 0.14.2 to rounding error, as long as no
#' ELBO change falls within rounding error of the tolerance.
#'
#' @noRd
NULL

#' @param X n x p matrix, no missing values
#' @param gram logical, work from the standardised Gram matrix
#' @param XtX crossprod(X), if already known
#' @noRd
.susie_prep <- function(X, gram = nrow(X) >= ncol(X) && .use_gram(ncol(X)),
                        XtX = NULL) {
    X <- as.matrix(X)
    storage.mode(X) <- "double"
    n <- nrow(X)
    p <- ncol(X)
    # susieR's compute_colstats with center = scale = TRUE; X reaches susie
    # as a sparse Matrix in univariate_susie, hence the sparse-matrix sd
    cm <- colMeans(X)
    csd_raw <- sqrt((colMeans(X^2) - cm^2) * n / (n - 1))
    csd <- csd_raw
    csd[csd == 0] <- 1
    d <- n * cm^2 + (n - 1) * csd_raw^2
    d <- (d - n * cm^2) / csd^2
    G <- NULL
    if (gram) {
        G <- (if (is.null(XtX)) crossprod(X) else XtX) - n * tcrossprod(cm)
        G <- G / csd
        G <- t(t(G) / csd)
    }
    list(X = X, n = n, p = p, cm = cm, csd = csd, d = d, G = G)
}

#' @param prep output of .susie_prep()
#' @param y numeric response
#' @return list(intercept, beta, converged, niter) with
#'   c(intercept, beta) equal to coef() of the susieR fit
#' @noRd
.susie_fit <- function(prep, y, scaled_prior_variance = 0.1,
                       max_iter = 100, tol = 1e-3) {
    X <- prep$X
    n <- prep$n
    p <- prep$p
    cm <- prep$cm
    csd <- prep$csd
    d <- prep$d
    G <- prep$G
    L <- min(10, p)
    y <- as.numeric(y)
    if (anyNA(y)) stop("Input y must not contain missing values.")

    lowerbound <- stats::var(y) / 1e4
    mean_y <- mean(y)
    y <- y - mean_y
    var_y <- stats::var(y)
    if (!(var_y > 0)) {
        stop("Residual variance sigma2 must be positive (is your var(Y) zero?).")
    }
    sigma2 <- var_y
    V <- scaled_prior_variance * var_y
    lpi <- log(rep(1 / p, p) + sqrt(.Machine$double.eps))

    # standardised products: Xs = t((t(X) - cm) / csd)
    Xty <- drop(crossprod(X, y)) / csd - cm / csd * sum(y)
    Xb <- function(b) drop(X %*% (b / csd)) - sum(cm * b / csd)

    alpha <- matrix(1 / p, L, p)
    mu <- matrix(0, L, p)
    mu2 <- matrix(0, L, p)
    KL <- rep(NA_real_, L)
    Eb <- matrix(0, L, p)                      # alpha * mu
    if (!is.null(G)) {
        Gb <- matrix(0, p, L)                  # Xs'Xs b_l
        GB <- numeric(p)                       # Xs'Xs sum_l b_l
    } else {
        Xbl <- matrix(0, n, L)                 # Xs b_l
        Xr <- numeric(n)
    }
    yty <- sum(y^2)

    prev_elbo <- -Inf
    converged <- FALSE
    for (iter in seq_len(max_iter)) {
        for (l in seq_len(L)) {
            if (!is.null(G)) {
                XtR <- Xty - (GB - Gb[, l])
            } else {
                R <- y - (Xr - Xbl[, l])
                XtR <- drop(crossprod(X, R)) / csd - cm / csd * sum(R)
            }
            betahat <- XtR / d
            shat2 <- sigma2 / d

            # single_effect_regression. lbf is dnorm(betahat, 0, sqrt(V +
            # shat2), log = TRUE) - dnorm(betahat, 0, sqrt(shat2), log = TRUE)
            # in closed form. A SNP with no variance (d = 0, shat2 = Inf)
            # gets log posterior odds 0, the prior variance and a mean of
            # (rounding error) 0.
            lbf <- -0.5 * log(1 + V / shat2) +
                0.5 * betahat^2 * V / (shat2 * (V + shat2))
            lpo <- lbf + lpi
            inf_idx <- is.infinite(shat2)
            lpo[inf_idx] <- 0
            mx <- max(lpo)
            w <- exp(lpo - mx)
            sw <- sum(w)
            a_l <- w / sw
            lbf_model <- log(sw) + mx
            post_var <- 1 / (1 / V + d / sigma2)
            post_mean <- post_var * XtR / sigma2

            alpha[l, ] <- a_l
            mu[l, ] <- post_mean
            mu2[l, ] <- post_var + post_mean^2
            eb <- a_l * post_mean
            Eb[l, ] <- eb
            # update_each_effect's KL, with the R'R and log 2 pi sigma2 terms
            # that cancel between its two halves removed
            KL[l] <- -lbf_model +
                (2 * sum(XtR * eb) - sum(d * a_l * mu2[l, ])) / (2 * sigma2)

            if (!is.null(G)) {
                gb <- drop(G %*% eb)
                GB <- GB - Gb[, l] + gb
                Gb[, l] <- gb
            } else {
                xb <- Xb(eb)
                Xr <- Xr - Xbl[, l] + xb
                Xbl[, l] <- xb
            }
        }

        # get_ER2
        B <- colSums(Eb)
        Eb2_d <- sum((alpha * mu2) %*% d)
        if (!is.null(G)) {
            rss <- yty - 2 * sum(Xty * B) + sum(B * GB)
            xb2 <- sum(Eb * t(Gb))
        } else {
            rss <- sum((y - Xr)^2)
            xb2 <- sum(Xbl^2)
        }
        ER2 <- rss + Eb2_d - xb2
        elbo <- -n / 2 * log(2 * pi * sigma2) - ER2 / (2 * sigma2) - sum(KL)
        if (is.infinite(elbo)) stop("get_objective() produced an infinite ELBO value")

        if (iter > 1) {
            diff <- elbo - prev_elbo
            # a decrease counts as converged too
            converged <- !is.na(diff) && diff < tol
        }
        prev_elbo <- elbo
        if (converged) break

        # estimate_residual_variance, bounded below as susie() does
        sigma2 <- max(lowerbound, ER2 / n)
    }
    if (!converged) {
        warning("IBSS algorithm did not converge in ", max_iter, " iterations!",
                call. = FALSE)
    }

    b <- colSums(alpha * mu) / csd
    list(intercept = mean_y - sum(cm * b),
         beta = b,
         converged = converged,
         niter = iter)
}

#' Coefficients (without intercept) of univariate_susie's fit for each column
#' of Y
#' @return p x ncol(Y) matrix
#' @noRd
.susie_coef <- function(X, Y, XtX = NULL) {
    prep <- .susie_prep(X, XtX = XtX)
    Y <- as.matrix(Y)
    vapply(seq_len(ncol(Y)), function(g) .susie_fit(prep, Y[, g])$beta,
           numeric(ncol(prep$X)))
}
