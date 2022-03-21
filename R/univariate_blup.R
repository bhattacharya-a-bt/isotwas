#' Univariate BLUP on all columns
#'
#' The function trains unvariate rrBLUP models individually for all
#' isoform transcripts
#'
#' @param X matrix, design matrix of SNP dosages
#' @param Y matrix, matrix of G isoform expression across columns
#' @param Omega matrix, precision matrix of Y
#' @param family character, glmnet glm family
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
#' @importFrom rrBLUP mixed.solve
#' @importFrom tibble tibble
#' @importFrom rlist list.append
#'
#' @export
univariate_blup <- function(X,
                            Y,
                            Omega,
                            scale = F,
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

    if (scale){
        Y = Y %*% (Omega %^% .5)
    }

    Y = as.data.frame(Y)
    if (!is.null(tx_names)){
        colnames(Y) = tx_names
    }
    set.seed(seed)
    train.folds = caret::createFolds(1:nrow(Y),
                                     k = nfolds,
                                     returnTrain = T)

    pred = matrix(ncol = ncol(Y),
                  nrow = nrow(Y))
    for (tr in 1:nfolds){

        Y.tr = as.matrix(Y[train.folds[[tr]],])
        X.tr = as.matrix(X[train.folds[[tr]],])
        X.test = as.matrix(X[-train.folds[[tr]],])
        blup.fit = apply(Y.tr,
                         FUN = rrBLUP::mixed.solve,
                         MARGIN = 2,
                         Z = X.tr,
                         K = t(X.tr) %*% X.tr / (ncol(X.tr) - 1))
        B.cur = sapply(blup.fit,
                       function(x) as.numeric(x$u))
        pred[-train.folds[[tr]],] = X.test %*% B.cur

    }

    r2.vec = sapply(1:ncol(Y),calc.r2,Y,pred)

    blup.fit = apply(Y,
                     FUN = rrBLUP::mixed.solve,
                     MARGIN = 2,
                     Z = X,
                     K = t(X) %*% X / (ncol(X) - 1))

    modelList = list()
    for (i in 1:length(blup.fit)){

        mod = tibble::tibble(SNP = colnames(X),
                             Weight = blup.fit[[i]]$u)
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
