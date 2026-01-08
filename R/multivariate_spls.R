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
                              seed){


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
  pred = matrix(ncol = ncol(Y),
                nrow = nrow(Y))

  set.seed(seed)
  max_K = min(ncol(Y),3)
  opt_model = spls::cv.spls(x = X,
                            y = Y,
                            fold = nfolds,
                            K = c(1:max_K),
                            eta = seq(.1,.9,.1),
                            scale.x=F,
                            scale.y=F,
                            plot.it = F)

  model_out = spls::spls(x = X,
                         y = Y,
                         K = opt_model$K.opt,
                         eta = opt_model$eta.opt,
                         scale.x=F,
                         scale.y=F)
  coef_mat = coef(model_out)

  # Define fold worker function
  run_fold <- function(tr) {
    Y.tr = Y[train.folds[[tr]],]
    X.tr = X[train.folds[[tr]],]
    test_idx = setdiff(1:nrow(Y), train.folds[[tr]])
    X.test = X[test_idx,]
    www = which(apply(X.tr,2,sd) != 0 &
                  apply(X.test,2,sd) != 0)
    X.tr = X.tr[,www]
    X.test = X.test[,www]

    fit <- spls::spls(x = X.tr,
                      y = Y.tr,
                      K = opt_model$K.opt,
                      eta = opt_model$eta.opt,
                      scale.x=F,
                      scale.y=F)

    list(test_idx = test_idx, preds = predict(fit, X.test))
  }

  # Run CV folds (parallel or sequential)
  if (par && !is.null(n.cores) && n.cores > 1) {
    fold_results <- parallel::mclapply(1:nfolds, run_fold, mc.cores = min(n.cores, nfolds))
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

  return(result)

}
