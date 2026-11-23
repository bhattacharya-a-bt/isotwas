#' Univariate BLUP on all columns (deprecated)
#'
#' Deprecated: use \code{\link{univariate_blup}}. \code{compute_blup} calls
#' \code{univariate_blup} and returns its result in the list format of
#' isotwas 2.0.0 and earlier. It will be removed in a future release.
#'
#' @inheritParams univariate_blup
#' @param seed int, random seed for drawing the folds
#'
#' @return list with one element per transcript, each a list of
#'   \code{Transcript} (name), \code{Model} (SNP weights), \code{R2},
#'   \code{P} and \code{Pred} (cross-validated predictions), as described
#'   for \code{\link{univariate_blup}}
#'
#' @export
compute_blup <- function(X,
                         Y,
                         Omega,
                         scale = F,
                         nfolds = 5,
                         verbose,
                         par = F,
                         n.cores = NULL,
                         tx_names = NULL,
                         seed){

  .Deprecated("univariate_blup")

  mod <- univariate_blup(X = X, Y = Y, Omega = Omega, scale = scale,
                         nfolds = nfolds, verbose = verbose, par = par,
                         n.cores = n.cores, tx_names = tx_names, seed = seed)

  lapply(unname(mod$transcripts), function(tm) {
    list(Transcript = tm$transcript_id,
         Model = tm$weights,
         R2 = tm$r2,
         P = tm$pvalue,
         Pred = tm$predicted)
  })
}
