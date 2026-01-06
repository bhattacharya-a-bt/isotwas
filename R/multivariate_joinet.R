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
#'
#' @return data frame of elastic net, lasso, and LMM based predictions
#'
#' @importFrom joinet cv.joinet
#' @importFrom joinet joinet
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
                                seed){

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

  if (par){
    doParallel::registerDoParallel(n.cores)
    set.seed(seed)
    models = joinet::cv.joinet(Y = Y,
                               X = X,
                               family="gaussian",
                               nfolds.ext=nfolds,
                               nfolds.int=10,
                               alpha.base=alpha,
                               alpha.meta=alpha,
                               parallel = T,
                               cvpred=T,
                               times=F,
                               trace.it = ifelse(verbose,1,0))
  } else  {
    set.seed(seed)
    models = joinet::cv.joinet(Y = Y,
                               X = X,
                               family="gaussian",
                               nfolds.ext=nfolds,
                               nfolds.int=10,
                               alpha.base=alpha,
                               alpha.meta=alpha,
                               parallel = F,
                               cvpred=T,
                               times=F,
                               trace.it = ifelse(verbose,1,0))
  }



  pred = models$cvpred
  r2.vec = unlist(sapply(1:ncol(Y),calc.r2,Y,pred)[1,])
  P = sapply(1:ncol(Y),function(x){
    stats::cor.test(Y[,x], pred[,x])$p.value
  })

  models = joinet::joinet(Y = Y,
                          X = X,
                          family="gaussian",
                          nfolds = 5,
                          alpha.base = alpha,
                          alpha.meta = alpha)

  # Build output using isotwas_model class
  transcripts <- list()
  for (i in 1:ncol(Y)){
    tx_name <- colnames(Y)[i]
    weights <- tibble::tibble(
      SNP = colnames(X),
      Weight = coef(models)$beta[, i]
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

  return(result)


}
