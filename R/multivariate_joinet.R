#' Multivariate stacked elastic net on all columns using joinet
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
#' @param nfolds.int int, folds for joinet's inner cross-validation.
#' @param folds optional cross-validation folds: a list with one vector of
#'   training-row indices per fold, or a vector of fold labels 1..K, one per
#'   sample. NULL (the default) draws them from \code{seed} as before;
#'   compute_isotwas passes every method the same folds.
#'   They are the external folds; the inner folds are drawn within them.
#'
#' @details
#' joinet is the most expensive method at moderate transcript counts, and most
#' of that time is its inner cross-validation. Lowering \code{nfolds.int} is
#' the first way to make it cheaper.
#'
#' @return data frame of elastic net, lasso, and LMM based predictions
#'
#' @importFrom pbapply pbapply
#' @importFrom tibble tibble
#' @importFrom rlist list.append
#' @importFrom stats cor.test
#'
#' @export
multivariate_joinet <- function(X,
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
                                nfolds.int = 10,
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

  # joinet::cv.joinet(Y, X, nfolds.ext = nfolds, nfolds.int = nfolds.int,
  # alpha.base = alpha, alpha.meta = alpha, cvpred = TRUE)$cvpred and then
  # joinet::joinet(Y, X, nfolds = 5, alpha.base = alpha, alpha.meta = alpha),
  # reproduced by .cv_joinet and .joinet_fit (joinet_fast.R) with the same
  # random draws, in the same order.
  set.seed(seed)
  Ym = as.matrix(Y)
  cores = if (par && !is.null(n.cores)) n.cores else 1L
  foldid.ext = if (is.null(folds)) NULL else .as_folds(folds, nrow(Ym))$foldid
  pred = .cv_joinet(X, Ym, nfolds.ext = nfolds, nfolds.int = nfolds.int,
                    alpha.base = alpha, alpha.meta = alpha, n.cores = cores,
                    foldid.ext = foldid.ext)
  r2.vec = unlist(sapply(1:ncol(Y),calc.r2,Y,pred)[1,])
  P = sapply(1:ncol(Y),function(x){
    stats::cor.test(Y[,x], pred[,x])$p.value
  })

  final = .joinet_fit(X, Ym, .joinet_folds(Ym[, 1], 5),
                      alpha.base = alpha, alpha.meta = alpha)
  beta = .joinet_coef(final)$beta

  # Build output using isotwas_model class
  transcripts <- list()
  for (i in 1:ncol(Y)){
    tx_name <- colnames(Y)[i]
    weights <- tibble::tibble(
      SNP = colnames(X),
      Weight = beta[, i]
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
    method = "joinet",
    transcripts = transcripts,
    n_samples = nrow(X),
    n_snps = ncol(X)
  )

  result <- .drop_null_fits(result, Y)

  return(result)


}
