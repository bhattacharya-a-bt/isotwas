#' Compute confirmation P-value for one gene
#'
#' The function runs a modified version of Shaffer's MSRB procedure
#' to control family-wise error rate for a given gene
#'
#' @param p numeric, vector of TWAS P-values for all isoforms of one gene
#' @param alpha numeric, significance level
#'
#' @return confirmation adjusted P-values for all isoforms of one gene
#'
#'
#' @export
p_confirm <- function(p,
                      alpha = 0.05){

    o <- order(p)
    n <- length(p)
    if (n==1) {
        adjustment=0
    } else {
        adjustment=c(n-1,(n-1):1)
    }
    pAdjusted <- p[o]*adjustment
    pAdjusted <- pmin(pAdjusted,1)
    pAdjusted <- cummax(pAdjusted)
    pBack <- vector(length=length(p))
    pBack[o] <- pAdjusted
    names(pBack) <- names(p)
    pBack[pBack > alpha] = 1
    return(pBack)

}


