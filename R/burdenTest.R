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
#' @importFrom stats pnorm
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

    if (is.null(Z)){
      if (any(is.null(beta,se))){
        stop('Please provide a column name for the Z-score or beta and SE.')
      }
      }

    colnames(sumStats)[which(colnames(sumStats) == chr)] = 'Chromosome'
    colnames(sumStats)[which(colnames(sumStats) == pos)] = 'Position'
    colnames(sumStats)[which(colnames(sumStats) == a1)] = 'A1'
    colnames(sumStats)[which(colnames(sumStats) == a2)] = 'A2'

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

      sumStats$SNP = paste(sumStats$Chromosome,sumStats$Position,sep=':')
      mod$SNP = paste(mod$Chromosome,mod$Position,sep=':')

    }

    tot = merge(mod,sumStats,by = 'SNP')

    if (nrow(tot) == 0){
      return('SNPs not found.')
    }

    tot$Z = ifelse(tot$A1.x == tot$A1.y,
                   tot$Z,
                   -1 * tot$Z)

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

    twasLD = as.numeric(tot$Weight %*% tot$Z) /
      sqrt(as.numeric(tot$Weight) %*% ld[tot$SNP,tot$SNP] %*% as.numeric(tot$Weight))
    twasLD = as.numeric(twasLD)
    P = 2*stats::pnorm(-abs(twasLD))

    if (P <= alpha){
      permutationLD = boot::boot(data = tot$Weight,
                                 statistic = calculateTWAS,
                                 R = nperms,
                                 sim = 'permutation',
                                 Z = tot$Z,
                                 LD = ld[tot$SNP,tot$SNP])
      permute.p = (nperms * mean(abs(permutationLD$t) >
                                   abs(permutationLD$t0)) + 1)/(nperms+1)
    } else {
      permute.p = 1}

  return(data.frame(Gene = gene,
                    Transcript = tot$Transcript[1],
                    Z = twasLD,
                    P = 2*pnorm(-abs(twasLD)),
                    permute.P = permute.p,
                    topSNP = tot$SNP[which.max(abs(tot$Z))],
                    topSNP.P = 2*pnorm(-abs(max(abs(tot$Z))))))
}
