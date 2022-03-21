#' Breiman and Friedman's curds and whey multivariate regression with CV
#'
#' The function trains the curds and whey multivariate regression with
#' cross-validated results
#'
#' @param X matrix, design matrix of SNP dosages
#' @param Y matrix, matrix of G isoform expression across columns
#' @param family character, glmnet glm family
#' @param scale logical, T/F to scale Y by Omega
#' @param alpha numeric, elastic net mixing parameter
#' @param nfolds int, number of CV folds
#' @param verbose logical
#' @param par logical, uses mclapply to parallelize model fit
#' @param n.cores int, number of parallel cores
#' @param tx_names vector, character vector of tx names - order of columns of Y
#' @param seed int, random seed
#'
#' @return data frame of elastic net, lasso, and LMM based predictions
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom tibble tibble
#' @importFrom rlist list.append
#' @importFrom caret createFolds
#'
#' @export
compute_curds_whey <- function(X,
                               Y,
                               family = 'gaussian',
                               alpha = 0.5,
                               nfolds = 5,
                               verbose,
                               par = F,
                               n.cores = NULL,
                               tx_names = NULL,
                               seed){

    if (!is.null(colnames(Y))){
        tx_names = colnames(Y)
    }

    Y = as.data.frame(Y)
    if (!is.null(tx_names)){
        colnames(Y) = tx_names
    }

    set.seed(seed)
    train.folds = caret::createFolds(1:nrow(Y),
                                     k = nfolds,
                                     returnTrain = T)

    pred = matrix(nrow = nrow(Y),
                  ncol = ncol(Y))
    for (tr in 1:nfolds){

        Y.tr = Y[train.folds[[tr]],]
        X.tr = X[train.folds[[tr]],]
        X.test = X[-train.folds[[tr]],]

        cur.B = curds_whey(X.tr,
                           Y.tr,
                           family = family,
                           alpha = alpha,
                           nfolds = nfolds,
                           verbose,
                           par = par,
                           n.cores = n.cores,
                           tx_names = tx_names)

        pred[-train.folds[[tr]],] = X.test %*% cur.B

    }
    pred[is.na(pred)] = 0
    r2.vec = sapply(1:ncol(Y),calc.r2,Y,pred)
    final.model = curds_whey(X,
                             Y,
                             family = family,
                             alpha = alpha,
                             nfolds = nfolds,
                             verbose,
                             par = par,
                             n.cores = n.cores,
                             tx_names = tx_names)


    modelList = list()
    for (i in 1:ncol(Y)){

        mod = tibble::tibble(SNP = colnames(X),
                             Weight = final.model[,i])
        mod = subset(mod,Weight != 0)
        modelList = rlist::list.append(modelList,
                                       list(Transcript = colnames(Y)[i],
                                            Model = mod,
                                            R2 = r2.vec[1,i],
                                            P = r2.vec[2,i],
                                            Pred = pred[,i]))

    }

    return(modelList)


}
