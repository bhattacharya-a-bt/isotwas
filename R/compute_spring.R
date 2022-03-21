#' Multivariate regression with covariance estimate using elastic net
#'
#' The function trains multivariate elastic-net with a given
#' covariance estimate
#'
#' @param X matrix, design matrix of SNP dosages
#' @param Y matrix, matrix of G isoform expression across columns
#' @param Omega matrix, precision matrix of Y
#' @param intercept logical, include intercept
#' @param normalize logical, normalize X and Y
#' @param nfolds int, number of CV folds
#' @param tol.in numeric, tolerance for objective difference
#' @param maxit.in int, maximum number of iterations
#' @param verbose logical
#' @param par logical, to parallelize
#' @param ncores int, number of cores
#' @param seed int, random seed
#'
#' @return CV MRCE fit
#'
#' @importFrom spring cv.spring
#' @importFrom MASS ginv
#' @importFrom tibble tibble
#' @importFrom rlist list.append
#' @importFrom caret createFolds
#'
#' @export
compute_spring = function(X,
                          Y,
                          Omega,
                          intercept = F,
                          normalize = T,
                          nfolds = 5,
                          tol.in,
                          maxit.in,
                          verbose = T,
                          par = F,
                          n.cores,
                          seed){

    if (!is.null(colnames(Y))){
        tx_names = colnames(Y)
    }

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

        www = which(apply(X.tr,2,sd) > 0 &
                        apply(X.test,2,sd) > 0)
        X.tr = X.tr[,www]
        X.test = X.test[,www]
        mod = suppressWarnings(spring::cv.spring(x = X.tr,
                                                 y = Y.tr,
                                                 cov = Omega,
                                                 intercept = intercept,
                                                 normalize = normalize,
                                                 threshold = tol.in,
                                                 max.iter = maxit.in,
                                                 verbose = ifelse(verbose,1,0),
                                                 mc.cores = n.cores))
        pred[-train.folds[[tr]],] =
            as.matrix(X.test %*% mod$cv.min$model$coef.regres)

    }

    r2.vec = sapply(1:ncol(Y),calc.r2,Y,pred)

    final.model = suppressWarnings(spring::cv.spring(x = X,
                                                     y = Y,
                                                     K = nfolds,
                                                     threshold = tol.in,
                                                     max.iter = maxit.in,
                                                     cov = MASS::ginv(Omega),
                                                     mc.cores = ifelse(par,
                                                                       n.cores,
                                                                       1)))
    Bhat = final.model$cv.min$model$coef.regres
    colnames(Y) = tx_names

    modelList = list()
    for (i in 1:ncol(Y)){

        mod = tibble::tibble(SNP = colnames(X),
                             Weight = as.numeric(Bhat[,i]))
        mod = subset(mod,Weight != 0)
        modelList = rlist::list.append(modelList,
                                       list(Transcript = colnames(Y)[i],
                                            Model = mod,
                                            R2 = r2.vec[i]))

    }

    return(modelList)


}
