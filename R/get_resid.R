#' Compute residuals for Z-scores
#'
#' The function runs a modified version of FOCUS's residual calculation
#'
#' @param zscores numeric, vector of Z-scores
#' @param swld numeric, LD matrix
#' @param wcor numeric, correlation matrix from estimate_cor
#'
#' @return residuals and residualized Z-scores
#'
#'
#' @export
get_resid = function(zscores, swld, wcor){
    m = nrow(wcor)
    p = ncol(swld)
    intercept = swld %*% rep(1,p)
    wcor_inv = MASS::ginv(as.matrix(wcor))
    rank = Matrix::rankMatrix(wcor_inv)[1]
    numer = t(intercept) %*% wcor_inv %*% zscores
    denom = t(intercept) %*% wcor_inv %*% intercept
    alpha = as.numeric(numer/denom)
    resid = zscores -
        (t(intercept) * alpha)
    s2 = (resid %*% wcor_inv %*%
              t(resid))/(rank-1)
    inter_se = sqrt(s2/denom)
    inter_z = alpha/inter_se
    return(list(resid,
                inter_z))
}
