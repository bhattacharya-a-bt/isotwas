#' Univariate elastic net on all columns
#'
#' The function trains unvariate elastic models individually for all
#' isoform transcripts
#'
#' @param X matrix, design matrix of SNP dosages
#' @param Y matrix, matrix of G isoform expression across columns
#' @param Omega matrix, precision matrix of Y
#' @param family character, glmnet glm family
#' @param scale logical, T/F to scale Y by Omega
#' @param alpha numeric, elastic net mixing parameter
#' @param nfolds int, number of CV folds
#' @param verbose logical
#' @param par logical, uses mclapply to parallelize model fit
#' @param n.cores int, number of parallel cores
#' @param tx_names vector, character vector of tx names in order of columns of Y
#' @param seed int, random seed
#' @param folds optional cross-validation folds: a list with one vector of
#'   training-row indices per fold, or a vector of fold labels 1..K, one per
#'   sample. NULL (the default) draws them from \code{seed} as before;
#'   compute_isotwas passes every method the same folds.
#'   Every transcript then uses these folds.
#'
#' @return data frame of elastic net, lasso, and LMM based predictions
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom parallel mclapply
#' @importFrom pbapply pbapply
#' @importFrom tibble tibble
#' @importFrom rlist list.append
#'
#' @export
univariate_elasticnet <- function(X,
                                  Y,
                                  Omega,
                                  family = 'gaussian',
                                  scale = F,
                                  alpha = 0.5,
                                  nfolds = 5,
                                  verbose,
                                  par = F,
                                  n.cores = NULL,
                                  tx_names = NULL,
                                  seed,
                                  folds = NULL){

    if (!is.null(colnames(Y))){
        tx_names = colnames(Y)
        }

    if (scale){
        Y = Y %*% (Omega %^% .5)
    }

    Y = as.data.frame(Y)
    if (!is.null(tx_names)){
        colnames(Y) = tx_names
        }

    foldid_given = NULL
    if (!is.null(folds)) {
        folds = .as_folds(folds, nrow(X))
        foldid_given = folds$foldid
        nfolds = folds$nfolds
    }

    if (family == 'gaussian'){
        # cv.glmnet(x = X, y = Y[, j], keep = TRUE, alpha = alpha) for each
        # transcript, reproduced by .cv_glmnet on compressed rows (see
        # glmnet_compress.R). Each call drew its own folds,
        # sample(rep(seq(nfolds), length = N)), in column order (after
        # set.seed(seed), or one seed per transcript with par = TRUE);
        # given folds (compute_isotwas's shared ones), every transcript uses
        # them. The full-data design serves every transcript; a fold design
        # serves its own transcript (so it is not kept in the run's memo),
        # or all of them with given folds.
        Ym = as.matrix(Y)
        N = nrow(X)
        n_fold = if (is.null(foldid_given)) N - ceiling(N / nfolds)
                 else N - max(tabulate(foldid_given))
        shared = NULL
        if (.glmnet_compress_on(n_fold, ncol(X))) {
            g = .gram_full_memo(X, Ym)
            shared = list(g = g, full = .gram_design(g))
            if (!is.null(foldid_given)) {
                shared$designs = .glmnet_cv_designs(X, Ym, foldid_given,
                                                    g = g, full = shared$full)
            }
        }
        fit_tx = function(j) {
            foldid = if (is.null(foldid_given)) {
                sample(rep(seq(nfolds), length = N))
            } else {
                foldid_given
            }
            designs = NULL
            if (!is.null(shared$designs)) {
                designs = shared$designs
            } else if (!is.null(shared)) {
                designs = .glmnet_cv_designs(X, Ym, foldid, g = shared$g,
                                             full = shared$full, keep = FALSE)
            }
            .cv_glmnet(X, Ym[, j], foldid, family = 'gaussian',
                       alpha = alpha, designs = designs, j = j)
        }
        if (par){
            models = .mclapply_seeded(seq_len(ncol(Ym)), FUN = fit_tx,
                                      mc.cores = n.cores, seed = seed)
        } else {
            set.seed(seed)
            models = if (verbose) {
                pbapply::pblapply(seq_len(ncol(Ym)), fit_tx)
            } else {
                lapply(seq_len(ncol(Ym)), fit_tx)
            }
        }
    } else {
        if (par){
            models = .mclapply_seeded(lapply(seq_len(ncol(Y)), function(i) Y[,i]),
                                      FUN = function(y) glmnet::cv.glmnet(
                                          x = X,
                                          y = y,
                                          nfolds = nfolds,
                                          foldid = foldid_given,
                                          keep = T,
                                          intercept = T,
                                          family = family,
                                          alpha = alpha),
                                      mc.cores = n.cores,
                                      seed = seed)
            } else if (verbose){
                set.seed(seed)
                models = pbapply::pbapply(Y,
                                          FUN = glmnet::cv.glmnet,
                                          MARGIN = 2,
                                          x = X,
                                          nfolds = nfolds,
                                          foldid = foldid_given,
                                          keep = T,
                                          intercept = T,
                                          family = family,
                                          alpha = alpha)
        } else {
            set.seed(seed)
            models = apply(Y,
                           FUN = glmnet::cv.glmnet,
                           2,
                           x = X,
                           nfolds = nfolds,
                           foldid = foldid_given,
                           keep = T,
                           intercept = T,
                           family = family,
                           alpha = alpha)
        }

    }

    # Build output using isotwas_model class
    transcripts <- list()
    for (i in 1:length(models)){
        tx_name <- colnames(Y)[i]
        weights <- tibble::tibble(
          SNP = colnames(X),
          Weight = coef(models[[i]], s = 'lambda.min')[-1]
        )
        weights <- subset(weights, Weight != 0)

        best.pred <- models[[i]]$fit.preval[, which.min(models[[i]]$cvm)]
        transcripts[[tx_name]] <- create_transcript_model(
          transcript_id = tx_name,
          weights = weights,
          r2 = .score_r2(as.numeric(Y[, i]), best.pred, pvalue = FALSE)$R2,
          pvalue = cor.test(as.numeric(Y[, i]), best.pred)$p.value,
          predicted = best.pred
        )
    }

    result <- create_isotwas_model(
      method = "univariate_enet",
      transcripts = transcripts,
      n_samples = nrow(X),
      n_snps = ncol(X)
    )

    result <- .drop_null_fits(result, Y)

    return(result)


}
