#' Group means and within-group residuals for the one-way design ~id
#'
#' Equal to the projection onto (or residual from) model.matrix(~id), without
#' forming the (R*N) x (R*N) hat matrix.
#'
#' @param Y.rep matrix, (R*N) x G isoform expression with replicates
#' @param id vector, length nrow(Y.rep), mapping each row to a sample
#' @param what character, 'fitted' for the per-sample means, 'residual' for
#'   the within-sample deviations
#'
#' @return matrix of the same dimension as Y.rep
#'
#' @noRd
.group_project <- function(Y.rep, id, what = c('fitted', 'residual')) {

    what = match.arg(what)

    if (nrow(Y.rep) != length(id)){
        stop("nrow(Y.rep) (", nrow(Y.rep), ") must equal length(id) (",
             length(id), ")")
    }

    grp = match(id, unique(id))
    means = rowsum(Y.rep, grp) / as.vector(table(grp))
    fitted = means[grp, , drop = FALSE]

    if (what == 'fitted') fitted else Y.rep - fitted
}

#' Compute precision matrix between isoforms
#'
#' The function computes the covariance and precision matrix between
#' columns of the isoform matrix
#'
#' @param Y matrix, matrix of G mean isoform expression across columns
#' @param Y.rep matrix, matrix of G isoform expression with replicates
#' @param R int, number of replicates
#' @param id vector, vector of sample ids showing rep to id
#' @param method character, one of 'replicates', 'residuals' or 'mean'; see
#'   Details.
#' @param verbose logical
#' @param nlambda int, number of matrices to produce in terms of sparsity
#'
#' @details
#' \code{'mean'} runs the graphical lasso on \code{Y}.
#'
#' \code{'replicates'} (the default) runs it on each sample's mean across
#' replicates in \code{Y.rep}, so within-sample variation is removed. With
#' \code{R = 1} it is identical to \code{'mean'}.
#'
#' \code{'residuals'} runs it on the within-sample deviations from those
#' means, the isoform-level technical covariance that replicates carry. It
#' needs at least two replicates per sample.
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
                          method = c('replicates','residuals','mean'),
                          nlambda = 15,
                          verbose = F){

    method = match.arg(method)

    if (!verbose){
        sink(file = tempfile())
        on.exit(sink(), add = TRUE)
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
        Y.fitted = .group_project(Y.rep, id, 'fitted')
        Omega_list = huge::huge(Y.fitted,
                                method = 'glasso',
                                nlambda=nlambda)
    } else if (method == 'residuals'){
        grp_sizes = as.vector(table(match(id, unique(id))))
        if (any(grp_sizes < 2)){
            stop("method = 'residuals' needs >= 2 replicates per sample; ",
                 sum(grp_sizes < 2), " of ", length(grp_sizes),
                 " samples have fewer, and their within-sample covariance is ",
                 "unidentifiable. Use method = 'mean' if replicates are ",
                 "genuinely unavailable.")
        }
        Y.resid = .group_project(Y.rep, id, 'residual')
        Omega_list = huge::huge(Y.resid,
                                method = 'glasso',
                                nlambda=nlambda)
    }

    return(Omega_list)


}
