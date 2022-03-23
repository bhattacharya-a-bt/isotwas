#' FOCUS-like gene-level finemapping
#'
#' The function runs probabilistic finemapping on the gene/isoform-level
#'
#' @param res.df data.frame, TWAS results that at least has Z scores and transcript names
#' @param z character, column name for Z scores
#' @param Omega matrix, matrix of SNP to transcript effect in same order as res.df
#' @param V matrix, LD matrix in the same order as Omega
#' @param max_enum numeric, max number of causals in the credible set
#' @param cutoff numeric, proportion of posterior explained by credible set
#'
#' @return data frame of res.df with pips and whether isoform is in credible set
#'
#' @export
twas_finemap <- function(res.df,
                         z = 'Z',
                         Omega,
                         V,
                         max_enum = 3,
                         cutoff = .9){

  colnames(res.df)[colnames(res.df) == z] = 'Z'
  zscores = res.df$Z
  m = length(zscores)
  wcor = estimate_cor(as.matrix(Omega),
                      as.matrix(V),
                      intercept=T)[[1]]
  diag(wcor) = 1
  wcor[is.na(wcor)] = 0
  swld = estimate_cor(as.matrix(Omega),
                      as.matrix(V),
                      intercept=T)[[2]]
  null_res = m * log(1 - 1e-3)
  marginal = m * log(1 - 1e-3)
  comb_list = list()
  for (n in 1:min(max_enum,length(zscores))){
    comb_list = c(comb_list,
                  combn(1:length(zscores),n,simplify=F))
  }
  pips = rep(0,length(zscores))
  zscores = get_resid(zscores,
                      as.matrix(swld),
                      as.matrix(wcor))[[1]]
  for (j in 1:length(comb_list)){
    subset = comb_list[[j]]
    local = bayes_factor(zscores,
                         idx_set = subset,
                         wcor = wcor)
    marginal = log(exp(local) + exp(marginal))
    for (idx in subset){
      if (pips[idx] == 0){
        pips[idx] = local
      } else {
        pips[idx] = log(exp(pips[idx]) + exp(local))
      }
    }
    print(pips)
    print(marginal)
  }

  pips = exp(pips - marginal)
  null_res = exp(null_res - marginal)
  res.df$pip = pips
  res.df = res.df[order(res.df$pip,decreasing = T),]
  npost = res.df$pip/sum(res.df$pip)
  csum = cumsum(npost)
  res.df$in_cred_set = csum < cutoff
  return(res.df)
}
