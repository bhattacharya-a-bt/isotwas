#' Multivariate regression with covariance estimate
#'
#' The function trains multivariate lasso with a given covariance estimate
#'
#' @param X matrix, design matrix of SNP dosages
#' @param Y matrix, matrix of G isoform expression across columns
#' @param lam2 numeric, lambda penalty for LASSO
#' @param weights matrix, weights for each predictor in X
#' @param Omega matrix, precision matrix of Y
#' @param tol.in numeric, tolerance for objective difference
#' @param maxit.in int, maximum number of iteractions
#' @param silent logical
#'
#' @return list of MRCE estimates
#'
#' @export
compute_fixed = function(X,
                         Y,
                         lam2,
                         Omega,
                         tol.in,
                         maxit.in,
                         silent){

    n=nrow(X)
    p=ncol(X)
    q=ncol(Y)
    if(!is.matrix(lam2)) {
        nlam=matrix(n*lam2, nrow=p, ncol=q)
    } else {
        nlam=n*lam2 }
    mx=apply(X, 2, mean)
    my=apply(Y, 2, mean)
    X=scale(X, center=mx, scale=FALSE)
    Y=scale(Y, center=my, scale=FALSE)
    xty=crossprod(X,Y)
    xtx=crossprod(X)
    tolmult=sum(crossprod(Y)*Omega)
    tol.in=tol.in*tolmult
    xtyom=xty%*%Omega
    new.B=rblasso(s=xtx, m=xtyom,
                  om=Omega, nlam=nlam,
                  n=n, B0=NULL,
                  soft=NULL, tol=tol.in,
                  maxit=maxit.in,
                  quiet=silent)$B
    muhat=as.numeric(my - crossprod(new.B, mx))
    return(list(Bhat=new.B, muhat=muhat,
                Omega=Omega, mx=mx, my=my))

}
