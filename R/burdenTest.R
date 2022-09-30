#' Compute weighted burden test
#'
#' The function takes in a gene expression model in MOSTWAS form
#' and GWAS summary statistics and
#' carries out the weighted burden Z-test for a trait
#'
#' @param mod data.frame, model for a given isoform
#' @param ld matrix, ld reference matrix
#' @param gene character, gene name
#' @param sumStats data frame, GWAS summary statistics
#' @param chr character, colnames in sumStats that keeps the chromosome
#' @param pos character, colnames in sumStats that keeps the position
#' @param a1 character, colnames in sumStats that keeps the ALT allele
#' @param a2 character, colnames in sumStats that keeps the REF allele
#' @param Z character, colnames in sumStats that keeps the Z score
#' @param beta character, colnames in sumStats that keeps the effect size
#' @param se character, colnames in sumStats that keeps the standard error
#' @param R2cutoff numeric, predictive R2 cutoff
#' @param alpha numeric, P-value threshold for permutation testing
#' @param nperms numeric, number of permutations
#' @param usePos logical, use SNP positions vs. SNP ids
#'
#' @return list of results for burden and permutation tests
#'
#' @importFrom boot boot
#'
#' @export
burdenTest <- function(mod,
                       ld,
                       gene,
                       sumStats,
                       chr,
                       pos,
                       a1,
                       a2,
                       Z = NULL,
                       beta = NULL,
                       se = NULL,
                       R2cutoff = 0.01,
                       alpha = 2.5e-6,
                       nperms = 1e3,
                       usePos = F){


  # if (all(is.null(c(Z,beta,se)))){
  #   stop('Please provide a column name for the Z-score or beta and SE.')
  # }
  #
  # if (is.null(Z) | any(is.null(c(beta,se)))){
  #   stop('Please provide a column name for the Z-score or beta and SE.')
  # }

  colnames(sumStats)[which(colnames(sumStats) == chr)] = 'Chromsome'
  colnames(sumStats)[which(colnames(sumStats) == pos)] = 'Position'

  if (!is.null(Z)){
    colnames(sumStats)[which(colnames(sumStats) == Z)] = 'Z'
  }

  if (!all(is.null(c(beta,se)))){
    colnames(sumStats)[which(colnames(sumStats) == beta)] = 'Beta'
    colnames(sumStats)[which(colnames(sumStats) == se)] = 'SE'
  }

  if (!'Z' %in% colnames(sumStats)){
    sumStats$Z = sumStats$Beta/sumStats$SE
  }

  if (mod$R2[1] <= R2cutoff){
    return(paste0('The isoform is not predicted at R2 > ',
                  R2cutoff))
  }

  if (usePos){

    sumStats$SNP = paste(sumStats$Chromsome,sumStats$Position,sep=':')
    mod$SNP = paste(mod$Chromosome,mod$Position,sep=':')

  }


  sumStats = subset(sumStats,SNP %in% mod$SNP)
  sumStats = sumStats[match(mod$SNP,
                            sumStats$SNP),]

  if (nrow(sumStats) == 0){
    return('SNPs not found.')
  }

  sumStats$Z = ifelse(sumStats$A1 == mod$A1,
                      sumStats$Z,
                      -1 * sumStats$Z)

  calculateTWAS <- function(effects,
                            Z,
                            LD,
                            indices){
    effects = effects[indices]
    twasZ = as.numeric(effects %*% Z)
    twasr2pred = as.numeric(effects %*% LD %*% effects)
    if (twasr2pred > 0){
      twas = as.numeric(twasZ/sqrt(twasr2pred))
    } else {
      twas = 0
    }
    return(twas)
  }

  twasLD = as.numeric(mod$Weight %*% sumStats$Z) /
    sqrt(as.numeric(mod$Weight) %*% ld[mod$SNP,mod$SNP] %*% as.numeric(mod$Weight))
  twasLD = as.numeric(twasLD)
  P = 2*pnorm(-abs(twasLD))

  if (P <= alpha){
    permutationLD = boot::boot(data = mod$Weight,
                               statistic = calculateTWAS,
                               R = nperms,
                               sim = 'permutation',
                               Z = sumStats$Z,
                               LD = ld[mod$SNP,mod$SNP])
    permute.p = (nperms * mean(abs(permutationLD$t) >
                                 abs(permutationLD$t0)) + 1)/(nperms+1)
  } else {
    permute.p = 1}

  return(data.frame(Gene = gene,
                    Transcript = mod$Transcript[1],
                    Z = twasLD,
                    P = 2*pnorm(-abs(twasLD)),
                    permute.P = permute.p,
                    topSNP = sumStats$SNP[which.max(abs(sumStats$Z))],
                    topSNP.P = 2*pnorm(-abs(max(abs(sumStats$Z))))))
}
