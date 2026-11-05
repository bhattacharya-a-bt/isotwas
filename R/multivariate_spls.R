#' Multivariate sparse partial least squares
#'
#' The function trains multivariate spls for all
#' isoform transcripts jointly
#'
#' @param X matrix, design matrix of SNP dosages
#' @param Y matrix, matrix of G isoform expression across columns
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
#' @return spls models
#'
#' @importFrom spls cv.spls
#' @importFrom spls spls
#' @importFrom tibble tibble
#' @importFrom rlist list.append
#' @importFrom caret createFolds
#' @importFrom stats sd
#'
#' @export
multivariate_spls <- function(X,
                              Y,
                              nfolds = 5,
                              verbose = T,
                              tx_names = NULL,
                              par,
                              n.cores,
                              seed,
                              folds = NULL){


  if (!is.null(colnames(Y))){
    tx_names = colnames(Y)
  }


  Y = as.data.frame(Y)
  if (!is.null(tx_names)){
    colnames(Y) = tx_names
  }

  X = X[,apply(X,2,stats::sd) != 0]

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

  # spls::cv.spls and spls::spls, reproduced exactly in .cv_spls and
  # .spls_fit (see spls_fast.R)
  Ym = as.matrix(Y)
  set.seed(seed)
  max_K = min(ncol(Y),3)
  opt_model = .cv_spls(x = X,
                       y = Ym,
                       fold = nfolds,
                       K = c(1:max_K),
                       eta = seq(.1,.9,.1))

  model_out = .spls_fit(.spls_prep(X, Ym),
                        K = opt_model$K.opt,
                        eta = opt_model$eta.opt)
  coef_mat = model_out$betahat

  # Define fold worker function
  run_fold <- function(tr) {
    Y.tr = Ym[train.folds[[tr]], , drop = FALSE]
    X.tr = X[train.folds[[tr]],]
    test_idx = setdiff(1:nrow(Y), train.folds[[tr]])
    X.test = X[test_idx,]
    www = which(apply(X.tr,2,sd) != 0 &
                  apply(X.test,2,sd) != 0)
    X.tr = X.tr[,www]
    X.test = X.test[,www]

    fit <- .spls_fit(.spls_prep(X.tr, Y.tr),
                     K = opt_model$K.opt,
                     eta = opt_model$eta.opt)

    # predict.spls: centre the selected columns of the new data
    A = fit$A
    newx = scale(X.test[, A, drop = FALSE], fit$meanx[A], rep(1, length(A)))
    preds = newx %*% fit$betahat[A, , drop = FALSE] +
      matrix(1, nrow(newx), 1) %*% fit$mu
    list(test_idx = test_idx, preds = preds)
  }

  # Run CV folds (parallel or sequential)
  if (par && !is.null(n.cores) && n.cores > 1) {
    fold_results <- .mclapply_seeded(1:nfolds, run_fold,
                                     mc.cores = min(n.cores, nfolds),
                                     seed = seed)
  } else {
    fold_results <- lapply(1:nfolds, run_fold)
  }

  # Aggregate predictions
  for (fr in fold_results) {
    pred[fr$test_idx,] <- fr$preds
  }
  r2.vec = unlist(sapply(1:ncol(Y),calc.r2,Y,pred)[1,])
  P = sapply(1:ncol(Y),function(x){
    cor.test(Y[,x], pred[,x])$p.value
  })

  # Build output using isotwas_model class
  transcripts <- list()
  for (i in 1:ncol(Y)){
    tx_name <- colnames(Y)[i]
    weights <- tibble::tibble(
      SNP = colnames(X),
      Weight = coef_mat[, i]
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
    method = "spls",
    transcripts = transcripts,
    n_samples = nrow(X),
    n_snps = ncol(X)
  )

  result <- .drop_null_fits(result, Y)

  return(result)

}
