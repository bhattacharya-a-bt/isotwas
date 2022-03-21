#' Compute correlation for fine-mapping
#'
#' The function runs a modified version of FOCUS's correlation matrix
#' estimation
#'
#' @param wmat numeric, matrix of weights
#' @param ldmat numeric, matrix of linkage disequilibrium
#' @param intercept logical, T/F if including an intercept
#'
#' @return correlation between GReX with/without intercept
#'
#'
#' @export
estimate_cor = function(wmat,ldmat,intercept=F){
    wcov = t(wmat) %*% ldmat %*% wmat
    scale = diag(1/sqrt(diag(wcov)))
    wcor = scale %*% wcov %*% scale
    if (intercept){
        inter = scale %*% t(wmat) %*% ldmat
        return(list(wcor,
                    inter))
    } else {
        return(list(wcor,NA))
    }
}
