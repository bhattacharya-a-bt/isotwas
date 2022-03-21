#' Multivariate regression with covariance estimate
#'
#' The function trains multivariate lasso with a given covariance estimate
#'
#' @param X matrix, design matrix of SNP dosages
#' @param Y matrix, matrix of G isoform expression across columns
#' @param lambda vector, lambda penalty vector for LASSO to tune on
#' @param Omega matrix, precision matrix of Y
#' @param nfolds int, number of CV folds
#' @param tol.in numeric, tolerance for objective difference
#' @param maxit.in int, maximum number of iteractions
#' @param verbose logical
#' @param seed int, random seed
#'
#' @return CV MRCE fit
#'
#' @export
compute_mrce = function(X,
                        Y,
                        lambda = NULL,
                        nlambda = 50,
                        Omega,
                        nfolds = 5,
                        tol.in,
                        maxit.in,
                        verbose,
                        seed){

    if (is.null(lambda)){
        lambda = 10^(seq(2,-30,length.out = nlambda))
        }

    set.seed(seed)
    train.folds = caret::createFolds(1:nrow(Y),
                                     k = nfolds,
                                     returnTrain = T)

    r2Mat = pMat = matrix(nrow = ncol(Y),
                   ncol = length(lambda))
    for (c in 1:length(lambda)){
        if (verbose) {print(c)}
        pred = matrix(nrow = nrow(Y),
                      ncol = ncol(Y))
        for (tr in 1:nfolds){
            Y.tr = Y[train.folds[[tr]],]
            X.tr = X[train.folds[[tr]],]
            X.test = X[-train.folds[[tr]],]

            www = which(apply(X.tr,2,sd) != 0 &
                            apply(X.test,2,sd) != 0)
            X.tr = X.tr[,www]
            X.test = X.test[,www]
            mod = compute_fixed(X = X.tr,
                                Y = Y.tr,
                                lam2 = lambda[c],
                                Omega = Omega,
                                tol.in = tol.in,
                                maxit.in = maxit.in,
                                silent = verbose)
            Bhat = mod$Bhat
            pred[-train.folds[[tr]],] = X.test %*% Bhat
        }
        pred[is.na(pred)] = 0

        r2 = as.numeric(as.matrix(sapply(1:ncol(pred),
                                         calc.r2,Y,pred)))[(1:(2*ncol(Y)))%%2 == 1]
        p = as.numeric(as.matrix(sapply(1:ncol(pred),
                                        calc.r2,Y,pred)))[(1:(2*ncol(Y)))%%2 == 0]
        r2[is.na(r2)] = 0
        p[is.na(p)] = 0
        r2Mat[,c] = r2
        pMat[,c] = p

    }

    final.model = compute_fixed(X = X,
                                Y = Y,
                                l = lambda[which.max(colMeans(r2Mat))],
                                Omega = Omega,
                                tol.in = tol.in,
                                maxit.in = maxit.in,
                                silent = verbose)


    modelList = list()
    for (i in 1:ncol(Y)){

        mod = tibble::tibble(SNP = colnames(X),
                             Weight = final.model$Bhat[,i])
        mod = subset(mod,Weight != 0)
        modelList = rlist::list.append(modelList,
                                       list(Transcript = colnames(Y)[i],
                                            Model = mod,
                                            R2 = r2Mat[i,
                                                       which.max(colMeans(r2Mat))],
                                            P = pMat[i,
                                                     which.max(colMeans(r2Mat))],
                                            Pred = pred[,i]))

    }

    return(modelList)




}
