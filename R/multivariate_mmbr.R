#' Multivariate SuSiE model
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
#' @param tx_names vector, character vector of tx names in order of columns of Y
#' @param L int, number of single effects
#' @param seed int, random seed
#'
#' @return data frame of mvSuSiE weights
#'
#' @importFrom mvsusieR mvsusie
#' @importFrom pbapply pbapply
#' @importFrom tibble tibble
#' @importFrom rlist list.append
#' @importFrom mvsusieR create_mash_prior
#' @importFrom stats cov
#'
#' @export
multivariate_mmbr <- function(X,
                              Y,
                              nfolds = 5,
                              verbose = T,
                              tx_names = NULL,
                              L = 10,
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
      print(tr)
        Y.tr = Y[train.folds[[tr]],]
        X.tr = X[train.folds[[tr]],]
        X.test = X[-train.folds[[tr]],]
        www = which(apply(X.tr,2,sd) != 0 &
                        apply(X.test,2,sd) != 0)
        X.tr = X.tr[,www]
        X.test = X.test[,www]
        prior_covar = mvsusieR::create_mash_prior(sample_data =
                                                  list(X=X.tr,
                                                       Y=Y.tr,
                                                       residual_variance =
                                                           stats::cov(Y.tr)),
                                              max_mixture_len=-1)
        m <- mvsusieR::mvsusie(X = X.tr,
                               Y = as.matrix(Y.tr),
                               prior_variance = prior_covar,
                               precompute_covariances = T,
                               L = L,
                               intercept = F,
                               standardize = T,
                               n_thread = ifelse(par,n.cores,1))

        pred[-train.folds[[tr]],] <- X.test %*% m$coef[-1,]
    }
    r2.vec = unlist(sapply(1:ncol(Y),calc.r2,Y,pred)[1,])
        P = sapply(1:ncol(Y),function(x){
          cor.test(Y[,x], pred[,x])$p.value
        })
        prior_covar = mvsusieR::create_mash_prior(sample_data =
                                                  list(X=X,
                                                       Y=Y,
                                                       residual_variance =
                                                           stats::cov(Y)),
                                              max_mixture_len=-1)
    m <- mvsusieR::mvsusie(X = X,
                           Y = as.matrix(Y),
                           prior_variance = prior_covar,
                           precompute_covariances = T,
                           L = L,
                           intercept = T,
                           standardize = T,
                           n_thread = ifelse(par,n.cores,1))

    modelList = list()
    for (i in 1:ncol(Y)){

        mod = tibble::tibble(SNP = colnames(X),
                             Weight = m$coef[-1,i])
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
