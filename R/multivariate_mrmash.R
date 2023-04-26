#' Multivariate MrMash model
#'
#' The function trains multivariate SuSiE for all
#' isoform transcripts jointly
#'
#' @param X matrix, design matrix of SNP dosages
#' @param Y matrix, matrix of G isoform expression across columns
#' @param nfolds int, number of CV folds
#' @param verbose logical
#' @param par logical, uses mclapply to parallelize model fit
#' @param n.cores int, number of parallel cores
#' @param tx_name vector, character vector of tx names in order of columns of Y
#'
#' @return mrmash models
#'
#' @importFrom mr.mash.alpha compute_univariate_sumstats
#' @importFrom mr.mash.alpha autoselect.mixsd
#' @importFrom mr.mash.alpha compute_canonical_covs
#' @importFrom mr.mash.alpha expand_covs
#' @importFrom mr.mash.alpha mr.mash
#' @importFrom pbapply pbapply
#' @importFrom tibble tibble
#' @importFrom rlist list.append
#' @importFrom caret createFolds
#'
#' @export
multivariate_mrmash <- function(X,
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

  X = X[,apply(X,2,sd) != 0]

  set.seed(seed)
  train.folds = caret::createFolds(1:nrow(Y),
                                   k = nfolds,
                                   returnTrain = T)
  pred = matrix(ncol = ncol(Y),
                nrow = nrow(Y))
  for (tr in 1:nfolds){
    Y.tr = Y[train.folds[[tr]],]
    X.tr = X[train.folds[[tr]],]
    X.test = X[-train.folds[[tr]],]
    www = which(apply(X.tr,2,sd) != 0 &
                  apply(X.test,2,sd) != 0)
    X.tr = X.tr[,www]
    X.test = X.test[,www]


    univ_sumstats <- mr.mash.alpha::compute_univariate_sumstats(X.tr, Y.tr,
                                                 standardize=TRUE,
                                                 standardize.response=FALSE,
                                                 mc.cores =
                                                   ifelse(par,n.cores,1))
    grid <- mr.mash.alpha::autoselect.mixsd(univ_sumstats,
                                            mult=sqrt(2))^2
    S0 <- mr.mash.alpha::compute_canonical_covs(ncol(Y.tr),
                                                singletons=TRUE,
                                 hetgrid=c(0, 0.25, 0.5, 0.75, 1))
    S0 <- mr.mash.alpha::expand_covs(S0,
                                     grid,
                                     zeromat=TRUE)

    fit <- mr.mash.alpha::mr.mash(X.tr,
                                  as.matrix(Y.tr),
                                  S0,
                                  update_V=TRUE,
                                  verbose = verbose,
                                  nthreads = ifelse(par,n.cores,1),
                                  version = 'Rcpp')

    pred[-train.folds[[tr]],] <- predict(fit,X.test)
  }
  r2.vec = unlist(sapply(1:ncol(Y),calc.r2,Y,pred)[1,])
  P = sapply(1:ncol(Y),function(x){
    cor.test(Y[,x], pred[,x])$p.value
  })


  univ_sumstats <-
    mr.mash.alpha::compute_univariate_sumstats(X,
                                               Y,
                                               standardize=TRUE,
                                               standardize.response=FALSE)
  grid <- mr.mash.alpha::autoselect.mixsd(univ_sumstats,
                                          mult=sqrt(2))^2
  S0 <- mr.mash.alpha::compute_canonical_covs(ncol(Y),
                                              singletons=TRUE,
                                              hetgrid=c(0, 0.25, 0.5, 0.75, 1))
  S0 <- mr.mash.alpha::expand_covs(S0,
                                   grid,
                                   zeromat=TRUE)

  fit <- mr.mash.alpha::mr.mash(X,
                                as.matrix(Y),
                                S0,
                                update_V=TRUE,
                                verbose = verbose,
                                nthreads = ifelse(par,n.cores,1),
                                version = 'Rcpp')

  modelList = list()
  for (i in 1:ncol(Y)){

    mod = tibble::tibble(SNP = colnames(X),
                         Weight = coef(fit)[-1,i])
    mod = subset(mod,Weight != 0)
    modelList = rlist::list.append(modelList,
                                   list(Transcript = colnames(Y)[i],
                                        Model = mod,
                                        R2 = r2.vec[i],
                                        P = P[i],
                                        Pred = pred[,i]))

  }


  return(modelList)

}
