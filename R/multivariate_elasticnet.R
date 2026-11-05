#' Multivariate elastic net on all columns with row-wise penalty
#'
#' The function trains multivariate elastic net models for all
#' isoform transcripts jointly
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
#'   All three alphas then use these folds.
#'
#' @return data frame of elastic net, lasso, and LMM based predictions
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom doParallel registerDoParallel
#' @importFrom pbapply pbapply
#' @importFrom tibble tibble
#' @importFrom rlist list.append
#'
#' @export
multivariate_elasticnet <- function(X,
                                    Y,
                                    Omega,
                                    scale = F,
                                    alpha = 0.5,
                                    nfolds = 5,
                                    verbose = T,
                                    par = F,
                                    n.cores = NULL,
                                    tx_names = NULL,
                                    seed,
                                    folds = NULL){

    if (!is.null(colnames(Y))){
        tx_names = colnames(Y)
    }

    if (scale){
        Y = Y %*% solve((Omega %^% .5))
    }

    Y = as.data.frame(Y)
    if (!is.null(tx_names)){
        colnames(Y) = tx_names
    }

    # cv.glmnet(x = X, y = Y, family = 'mgaussian', standardize = FALSE,
    # keep = TRUE) for alpha = 0, 0.5 and 1, reproduced by .cv_glmnet on
    # compressed rows (see glmnet_compress.R). Each cv.glmnet call drew its
    # own folds, sample(rep(seq(nfolds), length = N)), in this order. Given
    # folds (compute_isotwas's shared ones), all three alphas use them.
    set.seed(seed)
    Ym = as.matrix(Y)
    N = nrow(X)
    if (is.null(folds)) {
        foldids = lapply(1:3, function(i) sample(rep(seq(nfolds), length = N)))
    } else {
        folds = .as_folds(folds, N)
        foldids = rep(list(folds$foldid), 3)
        nfolds = folds$nfolds
    }
    # The full-data design serves all three alphas; a fold design serves
    # every alpha that uses its fold set.
    designs = list(NULL, NULL, NULL)
    if (.glmnet_compress_on(N - max(tabulate(foldids[[1]])), ncol(X),
                            ncol(Ym))) {
        g = .gram_full(X, Ym)
        full = .glmnet_design(.gram_subset(g))
        cv_designs = function(f) .glmnet_cv_designs(X, Ym, f, g = g,
                                                    full = full)
        designs = if (!is.null(folds)) {
            rep(list(cv_designs(foldids[[1]])), 3)
        } else {
            lapply(foldids, cv_designs)
        }
    }
    cores = if (par && !is.null(n.cores)) n.cores else 1L
    fit_alpha = function(i, a) {
        .cv_glmnet(X, Ym, foldids[[i]], family = 'mgaussian', alpha = a,
                   standardize = FALSE, designs = designs[[i]],
                   n.cores = cores)
    }
    models_0 = fit_alpha(1, 0)
    models_5 = fit_alpha(2, .5)
    models_1 = fit_alpha(3, 1)

    pred_5 = models_5$fit.preval[,,which.min(models_5$cvm)]
    pred_0 = models_0$fit.preval[,,which.min(models_0$cvm)]
    pred_1 = models_1$fit.preval[,,which.min(models_1$cvm)]

    r2.vec_0 = unlist(sapply(1:ncol(Y),calc.r2,Y,pred_0)[1,])
    r2.vec_5 = unlist(sapply(1:ncol(Y),calc.r2,Y,pred_5)[1,])
    r2.vec_1 = unlist(sapply(1:ncol(Y),calc.r2,Y,pred_1)[1,])

    r2_mat = matrix(nrow = 3,
                    ncol = ncol(Y))
    colnames(r2_mat) = colnames(Y)
    rownames(r2_mat) = c('alpha = 0','alpha = 0.5','alpha = 1')
    r2_mat[1,] = r2.vec_0
    r2_mat[2,] = r2.vec_5
    r2_mat[3,] = r2.vec_1

    P = sapply(1:ncol(Y),function(x){
      cor.test(Y[,x], pred_5[,x])$p.value
    })

    r2_best = vector('numeric',length = ncol(Y))
    best_alpha = vector('numeric',length = ncol(Y))
    for (t in 1:length(r2_best)){
      r2_best[t] = max(r2_mat[,t])
      best_alpha[t] = rownames(r2_mat)[which.max(r2_mat[,t])]
    }

    # Build output using isotwas_model class
    transcripts <- list()
    for (i in 1:ncol(Y)){
        tx_name <- colnames(Y)[i]

        if (best_alpha[i] == 'alpha = 0'){
          weights <- tibble::tibble(
            SNP = colnames(X),
            Weight = coef(models_0, s = 'lambda.min')[[i]][-1]
          )
          best_pred <- pred_0[, i]
        } else if (best_alpha[i] == 'alpha = 0.5'){
          weights <- tibble::tibble(
            SNP = colnames(X),
            Weight = coef(models_5, s = 'lambda.min')[[i]][-1]
          )
          best_pred <- pred_5[, i]
        } else {
          weights <- tibble::tibble(
            SNP = colnames(X),
            Weight = coef(models_1, s = 'lambda.min')[[i]][-1]
          )
          best_pred <- pred_1[, i]
        }

        weights <- subset(weights, Weight != 0)

        transcripts[[tx_name]] <- create_transcript_model(
          transcript_id = tx_name,
          weights = weights,
          r2 = r2_best[i],
          pvalue = P[i],
          predicted = best_pred
        )
    }

    result <- create_isotwas_model(
      method = "multi_enet",
      transcripts = transcripts,
      n_samples = nrow(X),
      n_snps = ncol(X)
    )

    result <- .drop_null_fits(result, Y)

    return(result)


}
