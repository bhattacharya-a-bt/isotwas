#' Univariate BLUP on all columns
#'
#' The function trains unvariate rrBLUP models individually for all
#' isoform transcripts
#'
#' @param X matrix, design matrix of SNP dosages
#' @param Y matrix, matrix of G isoform expression across columns
#' @param Omega matrix, precision matrix of Y
#' @param scale logical, T/F to scale Y by Omega
#' @param nfolds int, number of CV folds
#' @param verbose logical
#' @param par logical, uses mclapply to parallelize model fit
#' @param n.cores int, number of parallel cores
#' @param tx_names vector, character vector of tx names in order of columns of Y
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

    r2.vec = sapply(1:ncol(Y),function(x){
      summary((lm(Y[,x] ~ pred[,x])))$adj.r.sq
    })


    P = sapply(1:ncol(Y),function(x){
      cor.test(Y[,x], pred[,x])$p.value
    })

    blup.fit = apply(Y,
                     FUN = rrBLUP::mixed.solve,
                     MARGIN = 2,
                     Z = X,
                     K = t(X) %*% X / (ncol(X) - 1))

    # Build output using isotwas_model class
    transcripts <- list()
    for (i in 1:length(blup.fit)){
        tx_name <- colnames(Y)[i]
        weights <- tibble::tibble(
          SNP = colnames(X),
          Weight = blup.fit[[i]]$u
        )
        weights <- subset(weights, Weight != 0)

        transcripts[[tx_name]] <- create_transcript_model(
          transcript_id = tx_name,
          weights = weights,
          r2 = r2.vec[i],
          pvalue = P[i],
          predicted = pred[, i]
        )
    }

    result <- create_isotwas_model(
      method = "univariate_blup",
      transcripts = transcripts,
      n_samples = nrow(X),
      n_snps = ncol(X)
    )

    return(result)


}
