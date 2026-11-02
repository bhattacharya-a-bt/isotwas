# Deterministic fixtures shared across the test suite.
#
# Everything here depends on `stats` only, so the tests that exercise the
# numerical core (residuals, R2 scoring, model selection) run in environments
# where the heavier Bioconductor and modelling dependencies are unavailable.

#' Isoform expression with within-sample replicates
#'
#' @param N int, samples
#' @param G int, transcripts
#' @param R int, replicates per sample
#' @param sd numeric, within-sample replicate noise
#' @param unbalanced logical, if TRUE the group sizes are c(3, 2, 4) instead
#'   of R replicates for each of N samples
#' @return list(Y, Y.rep, id), with Y.rep stacked as rbind(rep_1, ..., rep_R)
#'   so that rows align with id
mk_reps <- function(N = 40, G = 6, R = 50, sd = 0.3, unbalanced = FALSE) {

    if (unbalanced) {
        id <- c(rep("A", 3), rep("B", 2), rep("C", 4))
        Y <- matrix(stats::rnorm(3 * G), 3, G,
                    dimnames = list(c("A", "B", "C"), paste0("tx", seq_len(G))))
        Y.rep <- matrix(stats::rnorm(length(id) * G), length(id), G,
                        dimnames = list(NULL, colnames(Y)))
        return(list(Y = Y, Y.rep = Y.rep, id = id))
    }

    tx <- paste0("tx", seq_len(G))
    samp <- paste0("S", seq_len(N))
    Y <- matrix(stats::rnorm(N * G), N, G, dimnames = list(samp, tx))
    reps <- lapply(seq_len(R), function(r) Y + matrix(stats::rnorm(N * G, sd = sd), N, G))
    Y.rep <- do.call(rbind, reps)
    colnames(Y.rep) <- tx

    list(Y = Y, Y.rep = Y.rep, id = rep(samp, times = R))
}

#' Genotype-like X and a Y generated from a known sparse coefficient matrix
#'
#' @param n int, samples
#' @param p int, SNPs
#' @param q int, transcripts
#' @param n_causal int, SNPs with nonzero effect
#' @param h2 numeric, proportion of variance explained by the genetic component
#' @return list(X, Y, B), all with dimnames set
mk_xy <- function(n = 100, p = 40, q = 3, n_causal = 5, h2 = 0.3) {

    maf <- stats::runif(p, 0.05, 0.5)
    X <- vapply(maf, function(f) stats::rbinom(n, 2, f), numeric(n))
    dimnames(X) <- list(paste0("S", seq_len(n)), paste0("rs", seq_len(p)))

    B <- matrix(0, p, q, dimnames = list(colnames(X), paste0("tx", seq_len(q))))
    causal <- sample(seq_len(p), n_causal)
    B[causal, ] <- stats::rnorm(n_causal * q)

    signal <- scale(X %*% B)
    noise <- matrix(stats::rnorm(n * q), n, q)
    Y <- sqrt(h2) * signal + sqrt(1 - h2) * noise
    dimnames(Y) <- list(rownames(X), colnames(B))

    list(X = X, Y = Y, B = B)
}

#' A plain list shaped like a fitted transcript model
#'
#' Carries every field create_isotwas_model() reads: weights, r2, pvalue and
#' the stored cross-validated predictions.
#'
#' @param r2 numeric, the cross-validated R2 to report
#' @param n_snps int, rows in the weight table; 0 produces an empty model
#' @param pvalue numeric, the reported p-value
#' @param n_pred int, length of the stored prediction vector
mk_tx <- function(r2, n_snps = 1L, pvalue = 0.05, n_pred = 10L) {
    weights <- if (n_snps == 0L) {
        data.frame(SNP = character(0), Weight = numeric(0))
    } else {
        data.frame(SNP = paste0("rs", seq_len(n_snps)),
                   Weight = rep(0.1, n_snps))
    }
    list(weights = weights,
         r2 = r2,
         pvalue = pvalue,
         predicted = stats::rnorm(n_pred))
}

#' A plain list shaped like an isotwas_model
#'
#' @param method character, the model's own $method string
#' @param r2s numeric, one cross-validated R2 per transcript
#' @param tx character, transcript names
#' @param n_snps int, recycled over transcripts; 0 produces an empty model
mk_model <- function(method, r2s, tx = paste0("t", seq_along(r2s)), n_snps = 1L) {
    n_snps <- rep_len(n_snps, length(r2s))
    list(method = method,
         transcripts = stats::setNames(
             lapply(seq_along(r2s), function(i) mk_tx(r2s[i], n_snps[i])), tx),
         n_samples = 100L,
         n_snps_total = 10L)
}

#' A plain list shaped like the result of a nested get_best() call
#'
#' Carries a $best_methods table, which is what distinguishes an ensemble from
#' a flat method when a label has to be resolved.
mk_ensemble <- function(r2s, leaves, tx = paste0("t", seq_along(r2s))) {
    mod <- mk_model("best", r2s, tx)
    mod$best_methods <- data.frame(transcript = tx,
                                   best_method = leaves,
                                   r2 = r2s,
                                   stringsAsFactors = FALSE)
    mod
}

#' The textbook projection residual, formed via an explicit hat matrix
#'
#' Deliberately the slow, obvious implementation: it is the reference the fast
#' path is checked against, so it must not share any code with it.
resid_hat <- function(Y.rep, id) {
    mm <- stats::model.matrix(~id)
    Y.rep - mm %*% solve(t(mm) %*% mm) %*% t(mm) %*% Y.rep
}

#' The same residual via QR, as a second independent reference
resid_qr <- function(Y.rep, id) {
    qr.resid(qr(stats::model.matrix(~id)), Y.rep)
}
