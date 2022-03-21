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
#' @param tx_names vector, character vector of tx names in order of columns of Y
#'
#' @return data frame of elastic net, lasso, and LMM based predictions
#'
#' @importFrom susieR susie
#' @importFrom estimatr lm_robust
#'
#' @export
finemap_regress <- function(X,
                            Y,
                            Y.rep,
                            R,
                            id,
                            verbose = F,
                            tx_names = NULL,
                            coverage = .9){

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

    susie.fit = list()
    for (i in 1:ncol(Y)){

        susie.fit =
            rlist::list.append(susie.fit,
                               susieR::susie(Matrix::Matrix(X,sparse = T),
                                             Y[,i],
                                             estimate_residual_variance = TRUE,
                                             estimate_prior_variance = FALSE,
                                             scaled_prior_variance = 0.1,
                                             verbose = verbose,
                                             coverage = coverage))

    }

    cs.all = sort(unique(unlist(sapply(susie.fit,
                                       function(x) unlist(x$sets$cs)))))
    P.fm = length(cs.all)
    if (P.fm <= 1){

        cs.all = c(cs.all,
                   sort(unique(unlist(sapply(susie.fit,
                                           function(x) which.max(x$pip))))))
        cs.all = c(cs.all,
                   sort(unique(unlist(sapply(susie.fit,
                                             function(x)
                                                 which(order(x$pip) == 2))))))
        cs.all = unique(cs.all)
        P.fm = length(cs.all)

    }

    X.design = do.call(rbind, replicate(R, X[,cs.all], simplify=FALSE))
    X.design = as.data.frame(X.design)
    total = as.data.frame(cbind(Y.rep,
                                X.design,
                                id))
    colnames(total)[1:G] = colnames(Y)
    fmla = as.formula(paste0('cbind(',
                             paste0('`',
                                    paste(colnames(Y),collapse = '`,`'),
                                    '`'),
                             ') ~ ',
                             paste0('`',paste(colnames(X[,cs.all]),
                                              collapse = '`+`'),'`')))

    regression = estimatr::lm_robust(formula = fmla,
                                     data = total,
                                     clusters = total$id)


    return(list(Coef = regression$coefficients,
                SE = regression$std.error,
                Error_Var = diag(cov(total[,1:G] - regression$fitted.values))))


}
