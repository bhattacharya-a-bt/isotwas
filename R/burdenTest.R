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
#' @param a1_mod character, colnames in model that keeps the ALT allele
#' @param a2_mod character, colnames in model that keeps the REF allele
#' @param snpName character, colnames in sumStats that keeps the SNP id
#' @param Z character, colnames in sumStats that keeps the Z score
#' @param beta character, colnames in sumStats that keeps the effect size
#' @param se character, colnames in sumStats that keeps the standard error
#' @param featureName character, colname in model that keeps the feature name
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
                       a1_mod = 'ALT',
                       a2_mod = 'REF',
                       snpName = 'SNP',
                       Z = NULL,
                       beta = NULL,
                       se = NULL,
                       featureName = 'Feature',
                       R2cutoff = 0.01,
                       alpha = 2.5e-6,
                       nperms = 1e3,
                       usePos = F){

    if (is.null(Z)){
      if (any(is.null(beta,se))){
        stop('Please provide a column name for the Z-score or beta and SE.')
      }
    }

   colnames(mod)[which(colnames(mod) == featureName)] = 'Feature'

   colnames(sumStats)[which(colnames(sumStats) == snpName)] = 'SNP'
   if (chr %in% colnames(sumStats)){
    colnames(sumStats)[which(colnames(sumStats) == chr)] = 'Chromosome'
   }
   if (pos %in% colnames(sumStats)){
    colnames(sumStats)[which(colnames(sumStats) == pos)] = 'Position'
   }
    colnames(sumStats)[which(colnames(sumStats) == a1)] = 'A1_GWAS'
    colnames(sumStats)[which(colnames(sumStats) == a2)] = 'A2_GWAS'
    colnames(sumStats)[which(colnames(mod) == a1_mod)] = 'A1_Mod'
    colnames(sumStats)[which(colnames(mod) == a2_mod)] = 'A2_Mod'

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

    tot$Z = ifelse(tot$A1_GWAS == tot$A1_Mod,
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
                    Feature = tot$Feature[1],
                    Z = twasLD,
                    P = 2*pnorm(-abs(twasLD)),
                    permute.P = permute.p,
                    topSNP = tot$SNP[which.max(abs(tot$Z))],
                    topSNP.P = 2*pnorm(-abs(max(abs(tot$Z))))))
}
