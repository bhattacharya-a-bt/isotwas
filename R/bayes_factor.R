#' Compute Bayes factors for TWAS Z-scores
#'
#' The function runs a modified version of FOCUS's Bayes factor
#' calculation
#'
#' @param zscores numeric, vector of Z-scores
#' @param idx_set character, ids of transcripts
#' @param wcor numeric, correlation matrix from estimate_cor
#' @param prior_chisq numeric, prior for chi-square mean
#' @param prb numeric, prior for causal configuration
#' @param use_log logical, use log of bayes factors
#'
#' @return bayes factor vector
#'
#'
#' @export
bayes_factor = function(zscores,
                        idx_set,
                        wcor,
                        prior_chisq = 40,
                        prb = 1e-3,
                        use_log = T){
    m = length(zscores)
    nc = length(idx_set)
    cur_chi2 = prior_chisq/nc
    cur_wcor = wcor[idx_set,idx_set]
    cur_zscors = zscores[idx_set]
    if (nc > 1){
        sss = svd(cur_wcor)
        cur_U = sss$u
        cur_EIG = sss$d
        rm(sss)
    } else {
        cur_U = 1
        cur_EIG = 1
    }
    scaled_chisq = cur_zscors^2

    cur_bf = .5 * -1*sum(log(1 + cur_chi2 %*% cur_EIG)) +
        .5 * sum((cur_chi2 / (1 + cur_chi2 %*% cur_EIG)) * scaled_chisq) +
        nc * log(prb) + (m-nc) * log(1-prb)
    if (use_log){
        return(cur_bf)
    } else {
        return(exp(cur_bf))
    }

}
