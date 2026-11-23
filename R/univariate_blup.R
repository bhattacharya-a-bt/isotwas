#' Univariate BLUP on all columns
#'
#' Fits a ridge-regression BLUP model to each isoform transcript (each column
#' of \code{Y}) separately, and scores it by cross-validation.
#'
#' Each transcript is fitted with the model
#' \eqn{y = 1 b + X u + e}{y = 1 b + X u + e}, where
#' \eqn{u \sim N(0, \sigma^2_u I)}{u ~ N(0, s2u I)},
#' \eqn{e \sim N(0, \sigma^2_e I)}{e ~ N(0, s2e I)} and the two variances
#' are estimated by REML, as
#' \code{rrBLUP::mixed.solve(y, Z = X, K = diag(ncol(X)))} would fit it.
#' Every SNP keeps a nonzero weight, except monomorphic ones.
#'
#' Models trained with isotwas 2.0.0 or earlier used a different prior, so
#' their weights differ; see NEWS for how to reproduce them.
#'
#' @param X matrix of SNP dosages, samples in rows and SNPs in columns, with
#'   SNP IDs as column names and no missing values
#' @param Y matrix of expression, samples in rows and one column per
#'   transcript
#' @param Omega matrix, precision matrix of \code{Y}; used only when
#'   \code{scale = TRUE}
#' @param scale logical, multiply \code{Y} by \code{Omega^(1/2)} before
#'   fitting
#' @param nfolds int, number of cross-validation folds, used when
#'   \code{folds} is \code{NULL}
#' @param verbose logical, not used by this method
#' @param par logical, not used by this method
#' @param n.cores int, not used by this method
#' @param tx_names character vector of transcript names, in the order of the
#'   columns of \code{Y}; ignored if \code{Y} has column names
#' @param seed int, random seed for drawing the folds when \code{folds} is
#'   \code{NULL}
#' @param folds optional cross-validation folds: a list with one vector of
#'   training-row indices per fold, or a vector of fold labels 1..K, one per
#'   sample. NULL (the default) draws them from \code{seed};
#'   compute_isotwas passes every method the same folds.
#'
#' @return an \code{isotwas_model} (see \code{\link{create_isotwas_model}})
#'   with, for each transcript, the SNP weights on the scale of \code{X}, the
#'   cross-validated R2 and p-value, and the cross-validated predictions
#'   (which include each training fold's intercept). The intercept of the
#'   final model is not stored. A transcript whose cross-validated prediction
#'   does not track its expression gets an empty set of weights and R2 0.
#'
#' @importFrom tibble tibble
#' @importFrom rlist list.append
#'
#' @export
univariate_blup <- function(X,
                            Y,
                            Omega,
                            scale = F,
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
    for (tr in 1:nfolds){

        Y.tr = as.matrix(Y[train.folds[[tr]], , drop = FALSE])
        X.tr = as.matrix(X[train.folds[[tr]], , drop = FALSE])
        X.test = as.matrix(X[-train.folds[[tr]], , drop = FALSE])
        fit.cur = .blup_fit(Y.tr, X.tr)
        pred[-train.folds[[tr]],] = sweep(X.test %*% fit.cur$weights, 2,
                                          fit.cur$intercept, "+")

    }

    r2.vec = sapply(1:ncol(Y),function(x){
      .score_r2(Y[,x], pred[,x], pvalue = FALSE)$R2
    })


    P = sapply(1:ncol(Y),function(x){
      cor.test(Y[,x], pred[,x])$p.value
    })

    U.full = .blup_fit(as.matrix(Y), as.matrix(X))$weights

    # Build output using isotwas_model class
    transcripts <- list()
    for (i in 1:ncol(U.full)){
        tx_name <- colnames(Y)[i]
        weights <- tibble::tibble(
          SNP = colnames(X),
          Weight = U.full[, i]
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
      method = "univariate_blup",
      transcripts = transcripts,
      n_samples = nrow(X),
      n_snps = ncol(X)
    )

    result <- .drop_null_fits(result, Y)

    return(result)


}
