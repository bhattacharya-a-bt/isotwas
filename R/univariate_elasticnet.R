#' Univariate elastic net on all columns
#'
#' The function trains unvariate elastic models individually for all
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
#' @importFrom glmnet cv.glmnet
#' @importFrom parallel mclapply
#' @importFrom pbapply pbapply
#' @importFrom tibble tibble
#' @importFrom rlist list.append
#'
#' @export
univariate_elasticnet <- function(X,
                                  Y,
                                  Omega,
                                  family = 'gaussian',
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

    if (par){
        set.seed(seed)
        models = parallel::mclapply(lapply(seq_len(ncol(Y)), function(i) x[,i]),
                                   FUN = glmnet::cv.glmnet,
                                   MARGIN = 2,
                                   x = X,
                                   nfolds = nfolds,
                                   keep = T,
                                   intercept = F,
                                   family = family,
                                   alpha = alpha,
                                   mc.cores = n.cores)
        } else if (verbose){
            set.seed(seed)
            models = pbapply::pbapply(Y,
                                      FUN = glmnet::cv.glmnet,
                                      MARGIN = 2,
                                      x = X,
                                      nfolds = nfolds,
                                      keep = T,
                                      intercept = F,
                                      family = family,
                                      alpha = alpha)
    } else {
        set.seed(seed)
        models = apply(Y,
                       FUN = glmnet::cv.glmnet,
                       2,
                       x = X,
                       nfolds = nfolds,
                       keep = T,
                       intercept = F,
                       family = family,
                       alpha = alpha)
    }

    modelList = list()
    for (i in 1:length(models)){

        mod = tibble::tibble(SNP = colnames(X),
                             Weight = coef(models[[i]],
                                           s = 'lambda.min')[-1])
        mod = subset(mod,Weight != 0)
        best.pred = models[[i]]$fit.preval[,which.min(models[[i]]$cvm)]
        reg = pred_r_squared(lm(Y[,i]~best.pred))
        modelList = rlist::list.append(modelList,
                                       list(Transcript = colnames(Y)[i],
                                            Model = mod,
                                            R2 = reg$r.squared,
                                            P = 1-pf(reg$fstatistic[1],
                                                     reg$fstatistic[2],
                                                     reg$fstatistic[3]),
                                            Pred = best.pred))

    }

    return(modelList)


}
