#' Univariate SuSiE on all columns
#'
#' The function trains unvariate SuSiE models individually for all
#' isoform transcripts
#'
#' @param X matrix, design matrix of SNP dosages
#' @param Y matrix, matrix of G isoform expression across columns
#' @param Omega matrix, precision matrix of Y
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
#'
#' @return data frame of elastic net, lasso, and LMM based predictions
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom parallel mclapply
#' @importFrom pbapply pbapply
#' @importFrom tibble tibble
#' @importFrom rlist list.append
#' @importFrom Matrix Matrix
#'
#' @export
univariate_susie <- function(X,
                             Y,
                             Omega,
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
    set.seed(seed)
    train.folds = caret::createFolds(1:nrow(Y),
                                     k = nfolds,
                                     returnTrain = T)
    if (!is.null(folds)) {
        folds = .as_folds(folds, nrow(Y))
        train.folds = folds$train
        nfolds = folds$nfolds
    }

    pred = matrix(ncol = ncol(Y),
                  nrow = nrow(Y))
    # X'X on all samples; each fold's is this minus its held-out rows
    X = as.matrix(X)
    XtX = if (nrow(X) >= ncol(X) && .use_gram(ncol(X))) crossprod(X) else NULL
    for (tr in 1:nfolds){

        Y.tr = as.matrix(Y[train.folds[[tr]],])
        X.tr = as.matrix(X[train.folds[[tr]],])
        X.test = as.matrix(X[-train.folds[[tr]],])
        # susieR::susie(X.tr, Y.tr[, i], estimate_residual_variance = TRUE,
        # estimate_prior_variance = FALSE, scaled_prior_variance = 0.1) for
        # every transcript, sharing the work on X.tr (see .susie_fit)
        B.cur = .susie_coef(X.tr, Y.tr,
                            XtX = if (!is.null(XtX)) XtX - crossprod(X.test))
        pred[-train.folds[[tr]],] = X.test %*% B.cur

    }

    #r2.vec = sapply(1:ncol(Y),calc.r2,Y,pred)
    r2.vec = sapply(1:ncol(Y),function(x){
      .score_r2(Y[,x], pred[,x], pvalue = FALSE)$R2
    })
    P = sapply(1:ncol(Y),function(x){
      cor.test(Y[,x], pred[,x])$p.value
    })

    B.full = .susie_coef(X, as.matrix(Y), XtX = XtX)

    # Build output using isotwas_model class
    transcripts <- list()
    for (i in 1:ncol(B.full)){
        tx_name <- colnames(Y)[i]
        weights <- tibble::tibble(
          SNP = colnames(X),
          Weight = B.full[, i]
        )
        weights <- subset(weights, Weight != 0)

        transcripts[[tx_name]] <- create_transcript_model(
          transcript_id = tx_name,
          weights = weights,
          r2 = r2.vec[i],
          pvalue = P[i],
          predicted = pred[, i]
        )
    }

    result <- create_isotwas_model(
      method = "univariate_susie",
      transcripts = transcripts,
      n_samples = nrow(X),
      n_snps = ncol(X)
    )

    result <- .drop_null_fits(result, Y)

    return(result)


}
