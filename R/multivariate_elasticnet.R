#' Multivariate elastic net on all columns with row-wise penalty
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
#' @importFrom glmnet cv.glmnet
#' @importFrom doParallel registerDoParallel
#' @importFrom pbapply pbapply
#' @importFrom tibble tibble
#' @importFrom rlist list.append
#'
#' @export
multivariate_elasticnet <- function(X,
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

    if (par){
        doParallel::registerDoParallel(n.cores)
        set.seed(seed)

        models_0 = glmnet::cv.glmnet(x = X,
                                     y = as.matrix(Y),
                                     nfolds = nfolds,
                                     parallel = par,
                                     family = 'mgaussian',
                                     keep = T,
                                     alpha = 0,
                                     intercept = T,
                                     standardize = F,
                                     trace.it = ifelse(verbose,1,0))

        models_5 = glmnet::cv.glmnet(x = X,
                                     y = as.matrix(Y),
                                     nfolds = nfolds,
                                     parallel = par,
                                     family = 'mgaussian',
                                     keep = T,
                                     alpha = .5,
                                     intercept = T,
                                     standardize = F,
                                     trace.it = ifelse(verbose,1,0))

        models_1 = glmnet::cv.glmnet(x = X,
                                   y = as.matrix(Y),
                                   nfolds = nfolds,
                                   parallel = par,
                                   family = 'mgaussian',
                                   keep = T,
                                   alpha = 1,
                                   intercept = T,
                                   standardize = F,
                                   trace.it = ifelse(verbose,1,0))
    } else  {
        set.seed(seed)

      models_0 = glmnet::cv.glmnet(x = X,
                                   y = as.matrix(Y),
                                   nfolds = nfolds,
                                   parallel = par,
                                   family = 'mgaussian',
                                   keep = T,
                                   alpha = 0,
                                   intercept = T,
                                   standardize = F,
                                   trace.it = ifelse(verbose,1,0))

      models_5 = glmnet::cv.glmnet(x = X,
                                   y = as.matrix(Y),
                                   nfolds = nfolds,
                                   parallel = par,
                                   family = 'mgaussian',
                                   keep = T,
                                   alpha = .5,
                                   intercept = T,
                                   standardize = F,
                                   trace.it = ifelse(verbose,1,0))

      models_1 = glmnet::cv.glmnet(x = X,
                                   y = as.matrix(Y),
                                   nfolds = nfolds,
                                   parallel = par,
                                   family = 'mgaussian',
                                   keep = T,
                                   alpha = 1,
                                   intercept = T,
                                   standardize = F,
                                   trace.it = ifelse(verbose,1,0))
    }

    pred_5 = models_5$fit.preval[,,which.min(models_5$cvm)]
    pred_0 = models_0$fit.preval[,,which.min(models_0$cvm)]
    pred_1 = models_1$fit.preval[,,which.min(models_1$cvm)]

    r2.vec_0 = unlist(sapply(1:ncol(Y),calc.r2,Y,pred_0)[1,])
    r2.vec_5 = unlist(sapply(1:ncol(Y),calc.r2,Y,pred_5)[1,])
    r2.vec_1 = unlist(sapply(1:ncol(Y),calc.r2,Y,pred_1)[1,])

    r2_mat = matrix(nrow = 3,
                    ncol = ncol(Y))
    colnames(r2_mat) = colnames(Y)
    rownames(r2_mat) = c('alpha = 0','alpha = 0.5','alpha = 1')
    r2_mat[1,] = r2.vec_0
    r2_mat[2,] = r2.vec_5
    r2_mat[3,] = r2.vec_1

    P = sapply(1:ncol(Y),function(x){
      cor.test(Y[,x], pred_5[,x])$p.value
    })

    r2_best = vector('numeric',length = ncol(Y))
    best_alpha = vector('numeric',length = ncol(Y))
    for (t in 1:length(r2_best)){
      r2_best[t] = max(r2_mat[,t])
      best_alpha[t] = rownames(r2_mat)[which.max(r2_mat[,t])]
    }

    modelList = list()
    for (i in 1:ncol(Y)){

        if (best_alpha[i] == 'alpha = 0'){

          mod = tibble::tibble(SNP = colnames(X),
                               Weight = coef(models_0,
                                             s = 'lambda.min')[[i]][-1])
          mod = subset(mod,Weight != 0)
          best.pred = models_0$fit.preval[,i,which.min(models_0$cvm)]
          modelList = rlist::list.append(modelList,
                                         list(Transcript = colnames(Y)[i],
                                              Model = mod,
                                              R2 = r2_best[i],
                                              P = P[i],
                                              Pred = pred_0[,i]))

        }

        if (best_alpha[i] == 'alpha = 0.5'){

          ccc = coef(models_5,
                     s = 'lambda.min')[[i]][-1]
          mod = tibble::tibble(SNP = colnames(X),
                               Weight = ccc)
          mod = subset(mod,Weight != 0)
          best.pred = models_5$fit.preval[,i,which.min(models_5$cvm)]
          modelList = rlist::list.append(modelList,
                                         list(Transcript = colnames(Y)[i],
                                              Model = mod,
                                              R2 = r2_best[i],
                                              P = P[i],
                                              Pred = pred_0[,i]))

        }

        if (best_alpha[i] == 'alpha = 1'){

          ccc = coef(models_1,
                     s = 'lambda.min')[[i]][-1]
          mod = tibble::tibble(SNP = colnames(X),
                               Weight = ccc)
          mod = subset(mod,Weight != 0)
          best.pred = models_1$fit.preval[,i,which.min(models_1$cvm)]
          modelList = rlist::list.append(modelList,
                                         list(Transcript = colnames(Y)[i],
                                              Model = mod,
                                              R2 = r2_best[i],
                                              P = P[i],
                                              Pred = pred_0[,i]))

        }


    }

    return(modelList)


}
