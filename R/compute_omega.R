#' Compute precision matrix between isoforms
#'
#' The function computes the covariance and precision matrix between
#' columns of the isoform matrix
#'
#' @param Y matrix, matrix of G mean isoform expression across columns
#' @param Y.rep matrix, matrix of G isoform expression with replicates
#' @param R int, number of replicates
#' @param id vector, vector of sample ids showing rep to id
#' @param method character, use Y.rep or Y
#' @param verbose logical
#'
#' @return covariance and/or precision matrix list in huge form
#'
#' @importFrom huge huge
#'
#' @export
compute_omega <- function(Y,
                          Y.rep,
                          R,
                          id,
                          method = c('replicates','mean'),
                          nlambda = 15,
                          verbose = F){

    if (!verbose){
        sink(file = tempfile())
    }
    G = ncol(Y)
    N = nrow(Y)

    if (is.null(R)){
        R = nrow(Y.rep)/nrow(Y)
    }

    if (method == 'mean'){
        Omega_list = huge::huge(Y,
                                method='glasso',
                                nlambda = nlambda)
    } else if (method == 'replicates'){
        ### RIGHT NOW THIS IS AN OLS WITH FIXED EFFECTS
        mm = model.matrix(~id)
        Y.resid = mm %*% solve(t(mm) %*% mm) %*% t(mm) %*% Y.rep
        Omega_list = huge::huge(Y.resid,
                                method = 'glasso',
                                nlambda=nlambda)
    }

    if (!verbose){
        sink()
    }

    return(Omega_list)


}
