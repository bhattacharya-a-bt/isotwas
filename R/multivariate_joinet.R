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
#' @param tx_name vector, character vector of tx names in order of columns of Y
#' @param seed int, random seed
#'
#' @return data frame of elastic net, lasso, and LMM based predictions
#'
#' @importFrom joinet cv.joinet
#' @importFrom joinet joinet
#' @importFrom pbapply pbapply
#' @importFrom tibble tibble
#' @importFrom rlist list.append
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

  models = joinet::cv.joinet(Y = Y,
                             X = X,
                             family="gaussian",
                             nfolds.ext=nfolds,
                             nfolds.int=10,
                             alpha.base=alpha,
                             alpha.meta=alpha,
                             cvpred=T,
                             times=F)

  pred = models$cvpred
  r2.vec = unlist(sapply(1:ncol(Y),calc.r2,Y,pred)[1,])
  P = sapply(1:ncol(Y),function(x){
    cor.test(Y[,x], pred[,x])$p.value
  })

  models = joinet::joinet(Y = Y,
                          X = X,
                          family="gaussian",
                          nfolds = 5,
                          alpha.base = alpha,
                          alpha.meta = alpha)

  modelList = list()
  for (i in 1:ncol(Y)){

    mod = tibble::tibble(SNP = colnames(X),
                         Weight = coef(models)$beta[,i])
    mod = subset(mod,Weight != 0)
    best.pred = pred[,i]
    reg = summary(lm(Y[,i]~best.pred))
    modelList = rlist::list.append(modelList,
                                   list(Transcript = colnames(Y)[i],
                                        Model = mod,
                                        R2 = r2.vec[i],
                                        P = P[i],
                                        Pred = pred[,i]))

  }

  return(modelList)


}
