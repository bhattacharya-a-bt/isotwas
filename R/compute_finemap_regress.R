#' Fine-mapping and clustering SEs with SUR
#'
#' The function feature selects with SuSiE and runs SUR with cluster-robust
#' SEs
#'
#' @param X matrix, design matrix of SNP dosages
#' @param Y matrix, matrix of G mean isoform expression across columns
#' @param Y.rep matrix, matrix of G isoform expression with replicates
#' @param R int, number of replicates
#' @param id vector, vector of sample ids showing rep to id
#' @param verbose logical
#' @param tx_names vector, character vector of tx names - order of columns of Y
#' @param seed int, random seed
#'
#' @return isoTWAS model list
#'
#' @importFrom susieR susie
#' @importFrom caret createFolds
#' @importFrom tibble tibble
#' @importFrom rlist list.append
#'
#' @export
compute_finemap_regress <- function(X,
                                    Y,
                                    Y.rep,
                                    R,
                                    id,
                                    nfolds,
                                    verbose = F,
                                    tx_names = NULL,
                                    coverage = .9,
                                    seed){

    G = ncol(Y)
    N = nrow(Y)
    P = ncol(X)
    if (is.null(R)){
        R = mean(table(as.factor(id)))
    }

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

    pred = matrix(ncol = ncol(Y),
                  nrow = nrow(Y))
    for (tr in 1:nfolds){

        Y.tr = Y[train.folds[[tr]],]
        X.tr = X[train.folds[[tr]],]
        X.test = X[-train.folds[[tr]],]
        colnames(X.test) = sapply(colnames(X.test),
                                  function(x) paste0('`',x,'`'))
        id.tr = as.character(id[id %in% levels(id)[train.folds[[tr]]]])
        id.tr = as.factor(id.tr)
        Y.rep.tr = Y.rep[id %in% levels(id)[train.folds[[tr]]],]
        B.cur = finemap_regress(X = X.tr,
                                Y = Y.tr,
                                Y.rep = Y.rep.tr,
                                R = R,
                                id = id.tr,
                                verbose = F,
                                tx_names = NULL,
                                coverage = coverage)$Coef[-1,]
        colnames(X.test) = sapply(strsplit(colnames(X.test),'`'),
                                  function(x) x[2])
        pred[-train.folds[[tr]],] = X.test[,rownames(B.cur)] %*% B.cur

    }
    pred[is.na(pred)] = 0

    r2.vec = sapply(1:ncol(Y),calc.r2,Y,pred)
    final.model = finemap_regress(X = X,
                                  Y = Y,
                                  Y.rep = Y.rep,
                                  R = R,
                                  id = id,
                                  verbose = F,
                                  tx_names = NULL,
                                  coverage = coverage)

    modelList = list()
    for (i in 1:ncol(Y)){

        mod = tibble::tibble(SNP = rownames(final.model$Coef),
                             Weight = as.numeric(final.model$Coef[,i]),
                             SE = as.numeric(final.model$SE[,i]))
        mod = subset(mod,Weight != 0)
        modelList = rlist::list.append(modelList,
                                       list(Transcript = colnames(Y)[i],
                                            Model = mod,
                                            Error_Var = final.model$Error_Var,
                                            R2 = c(r2.vec[1,i]),
                                            P = c(r2.vec[2,i]),
                                            Pred = pred[,i]))

    }

    return(modelList)


}
