#' Breiman and Friedman's curds and whey multivariate regression
#'
#' The function trains the curds and whey multivariate regression
#'
#' @param X matrix, design matrix of SNP dosages
#' @param Y matrix, matrix of G isoform expression across columns
#' @param family character, glmnet glm family
#' @param alpha numeric, elastic net mixing parameter
#' @param nfolds int, number of CV folds
#' @param verbose logical
#' @param par logical, uses mclapply to parallelize model fit
#' @param n.cores int, number of parallel cores
#' @param tx_names vector, character vector of tx names in order of columns of Y
#'
#' @return data frame of elastic net, lasso, and LMM based predictions
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom parallel mclapply
#' @importFrom pbapply pbapply
#' @importFrom tibble tibble
#' @importFrom rlist list.append
#' @importFrom stats cancor
#' @importFrom stats runif
#'
#' @keywords internal
#' @export
curds_whey <- function(X,
                       Y,
                       family = 'gaussian',
                       alpha = 0.5,
                       nfolds = 5,
                       verbose,
                       par = F,
                       n.cores = NULL,
                       tx_names = NULL){

    if (!is.null(colnames(Y))){
        tx_names = colnames(Y)
    }

    Y = as.data.frame(Y)
    if (!is.null(tx_names)){
        colnames(Y) = tx_names
    }

    Y = as.matrix(scale(Y))
    #X = as.matrix(scale(X))
    cancor.all <- stats::cancor(X, Y)
    cancor.cor <- cancor.all$cor
    cancor.y <- cancor.all$ycoef
    yPrime <- as.matrix(Y) %*% cancor.y
    yPrime.mod = apply(yPrime,
                       FUN = glmnet::cv.glmnet,
                       2,
                       x = X,
                       nfolds = nfolds,
                       alpha = alpha,
                       family = family)
    B.yPrime.mod = sapply(yPrime.mod,
                          function(x) as.numeric(coef(x,s='lambda.min'))[-1])
    yhat.Prime = X %*% B.yPrime.mod
    r = ncol(X) / nrow(X)
    di = rep(0, ncol(Y))
    di = {(1 - r) * ({cancor.cor^2} - r)} /
        + {({1 - r}^2) * (cancor.cor^2) + ({r^2} * {1 - cancor.cor^2}) }
    for(i in 1:ncol(Y)){
        di[i] <- max(di[i], 0)
    }
    D <- diag(di)
    yhatstar <- yhat.Prime %*% D
    yFinal <- yhatstar %*% solve(cancor.y)
    www = which(apply(yFinal,2,sd) == 0)
    yFinal[,www] = stats::runif(length(www)*nrow(yFinal),
                         -0.1,.1)

    models = apply(yFinal,
                   FUN = glmnet::cv.glmnet,
                   2,
                   x = X,
                   nfolds = nfolds,
                   alpha = alpha,
                   family = family)


    return(sapply(models,
                  function(x) as.numeric(coef(x,s='lambda.min'))[-1]))


}
