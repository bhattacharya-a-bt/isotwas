#' Compute weighted burden test
#'
#' The function takes in a gene expression model (either legacy format or
#' new isotwas_model format) and GWAS summary statistics and carries out
#' the weighted burden Z-test for a trait.
#'
#' @param mod Model for a given isoform. Can be:
#'   - A data.frame in legacy format with columns: Feature, SNP, Weight, R2, etc.
#'   - A single transcript model from isotwas_model$transcripts$<name>
#'   - An isotwas_model object (will use first transcript)
#' @param ld matrix, ld reference matrix with SNP names as row/colnames
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
#' @return data.frame with results for burden and permutation tests
#'
#' @importFrom boot boot
#' @importFrom stats pnorm
#'
#' @export
burdenTest <- function(mod,
                       ld,
                       gene,
                       sumStats,
                       chr = NULL,
                       pos = NULL,
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
                       usePos = FALSE){

  if (is.null(Z)){
    if (any(is.null(c(beta,se)))){
      stop('Please provide a column name for the Z-score or beta and SE.')
    }
  }

  # Convert new model format to legacy format if needed
  mod <- .convert_model_for_burden(mod, featureName)

  # Handle case where conversion failed

  if (is.character(mod) && length(mod) == 1) {
    return(mod)  # Return error message
  }

  colnames(mod)[which(colnames(mod) == featureName)] = 'Feature'

  colnames(sumStats)[which(colnames(sumStats) == snpName)] = 'SNP'
  if (!is.null(chr) && chr %in% colnames(sumStats)){
    colnames(sumStats)[which(colnames(sumStats) == chr)] = 'Chromosome'
  }
  if (!is.null(pos) && pos %in% colnames(sumStats)){
    colnames(sumStats)[which(colnames(sumStats) == pos)] = 'Position'
  }
  colnames(sumStats)[which(colnames(sumStats) == a1)] = 'A1_GWAS'
  colnames(sumStats)[which(colnames(sumStats) == a2)] = 'A2_GWAS'

  if (a1_mod %in% colnames(mod)) {
    colnames(mod)[which(colnames(mod) == a1_mod)] = 'A1_Mod'
  }
  if (a2_mod %in% colnames(mod)) {
    colnames(mod)[which(colnames(mod) == a2_mod)] = 'A2_Mod'
  }

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
    sumStats$SNP = paste(sumStats$Chromosome, sumStats$Position, sep=':')
    mod$SNP = paste(mod$Chromosome, mod$Position, sep=':')
  }

  tot = merge(mod, sumStats, by = 'SNP')

  # making sure all alleles are in the same case
  tot$A1_GWAS = toupper(tot$A1_GWAS)
  tot$A2_GWAS = toupper(tot$A2_GWAS)

  if ('A1_Mod' %in% colnames(tot)) {
    tot$A1_Mod = toupper(tot$A1_Mod)
    tot$A2_Mod = toupper(tot$A2_Mod)

    # Flip Z-score if alleles are swapped
    tot$Z = ifelse(tot$A1_GWAS == tot$A1_Mod,
                   tot$Z,
                   -1 * tot$Z)
  }

  if (nrow(tot) == 0){
    return('SNPs not found.')
  }

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

  # Check that SNPs are in LD matrix
  snps_in_ld <- tot$SNP[tot$SNP %in% rownames(ld)]
  if (length(snps_in_ld) == 0) {
    return('No SNPs found in LD matrix.')
  }
  tot <- tot[tot$SNP %in% snps_in_ld, ]

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
    permute.p = 1
  }

  return(data.frame(Gene = gene,
                    Feature = tot$Feature[1],
                    Z = twasLD,
                    P = 2*pnorm(-abs(twasLD)),
                    permute.P = permute.p,
                    topSNP = tot$SNP[which.max(abs(tot$Z))],
                    topSNP.P = 2*pnorm(-abs(max(abs(tot$Z))))))
}


#' Convert model to legacy format for burden test
#'
#' Internal function to convert isotwas_model or transcript model to
#' the data.frame format expected by burdenTest
#'
#' @param mod Model object (various formats supported)
#' @param featureName Column name for feature/transcript ID
#'
#' @return data.frame in legacy format
#'
#' @keywords internal
.convert_model_for_burden <- function(mod, featureName = 'Feature') {


  # If it's already a data.frame with required columns, return as-is
  if (is.data.frame(mod) && 'SNP' %in% colnames(mod) && 'Weight' %in% colnames(mod)) {
    return(mod)
  }

  # If it's an isotwas_model, extract first transcript
  if (inherits(mod, 'isotwas_model')) {
    if (length(mod$transcripts) == 0) {
      return('Model has no transcripts.')
    }
    tx_name <- names(mod$transcripts)[1]
    mod <- mod$transcripts[[tx_name]]
  }

  # If it's a single transcript model from isotwas_model$transcripts$<name>
  if (is.list(mod) && all(c('weights', 'r2', 'transcript_id') %in% names(mod))) {
    weights_df <- mod$weights

    if (nrow(weights_df) == 0) {
      return('Model has no non-zero weights.')
    }

    # Build data.frame in legacy format
    result <- data.frame(
      Feature = mod$transcript_id,
      SNP = weights_df$SNP,
      Weight = weights_df$Weight,
      R2 = mod$r2,
      stringsAsFactors = FALSE
    )

    # Add allele columns if present in weights
    if ('ALT' %in% colnames(weights_df)) {
      result$ALT <- weights_df$ALT
    }
    if ('REF' %in% colnames(weights_df)) {
      result$REF <- weights_df$REF
    }

    return(result)
  }

  # If we get here, format not recognized
  return('Unrecognized model format.')
}


#' Run burden test on all transcripts in an isotwas_model
#'
#' Convenience function to run burdenTest on all transcripts in an
#' isotwas_model or isotwas_result object.
#'
#' @param model An isotwas_model or isotwas_result object
#' @param ld matrix, ld reference matrix with SNP names as row/colnames
#' @param gene character, gene name
#' @param sumStats data frame, GWAS summary statistics
#' @param chr character, colnames in sumStats that keeps the chromosome
#' @param pos character, colnames in sumStats that keeps the position
#' @param a1 character, colnames in sumStats that keeps the ALT allele
#' @param a2 character, colnames in sumStats that keeps the REF allele
#' @param snpName character, colnames in sumStats that keeps the SNP id
#' @param Z character, colnames in sumStats that keeps the Z score
#' @param beta character, colnames in sumStats that keeps the effect size
#' @param se character, colnames in sumStats that keeps the standard error
#' @param R2cutoff numeric, predictive R2 cutoff
#' @param alpha numeric, P-value threshold for permutation testing
#' @param nperms numeric, number of permutations
#' @param usePos logical, use SNP positions vs. SNP ids
#' @param verbose logical, print progress
#'
#' @return data.frame with burden test results for all transcripts
#'
#' @export
burdenTest_all <- function(model,
                           ld,
                           gene,
                           sumStats,
                           chr = NULL,
                           pos = NULL,
                           a1,
                           a2,
                           snpName = 'SNP',
                           Z = NULL,
                           beta = NULL,
                           se = NULL,
                           R2cutoff = 0.01,
                           alpha = 2.5e-6,
                           nperms = 1e3,
                           usePos = FALSE,
                           verbose = FALSE) {

  # Extract transcripts from model
  if (inherits(model, 'isotwas_result')) {
    transcripts <- model$best_models$transcripts
  } else if (inherits(model, 'isotwas_model')) {
    transcripts <- model$transcripts
  } else {
    stop('model must be an isotwas_model or isotwas_result object')
  }

  tx_names <- names(transcripts)
  results_list <- list()

  for (i in seq_along(tx_names)) {
    tx_name <- tx_names[i]
    if (verbose) {
      cat(sprintf("Testing %s (%d/%d)\n", tx_name, i, length(tx_names)))
    }

    tx_mod <- transcripts[[tx_name]]

    result <- tryCatch({
      burdenTest(
        mod = tx_mod,
        ld = ld,
        gene = gene,
        sumStats = sumStats,
        chr = chr,
        pos = pos,
        a1 = a1,
        a2 = a2,
        snpName = snpName,
        Z = Z,
        beta = beta,
        se = se,
        R2cutoff = R2cutoff,
        alpha = alpha,
        nperms = nperms,
        usePos = usePos
      )
    }, error = function(e) {
      data.frame(
        Gene = gene,
        Feature = tx_name,
        Z = NA,
        P = NA,
        permute.P = NA,
        topSNP = NA,
        topSNP.P = NA,
        Error = e$message,
        stringsAsFactors = FALSE
      )
    })

    # Handle character returns (error messages)
    if (is.character(result)) {
      result <- data.frame(
        Gene = gene,
        Feature = tx_name,
        Z = NA,
        P = NA,
        permute.P = NA,
        topSNP = NA,
        topSNP.P = NA,
        Note = result,
        stringsAsFactors = FALSE
      )
    }

    results_list[[tx_name]] <- result
  }

  # Combine all results
  results_df <- do.call(rbind, results_list)
  rownames(results_df) <- NULL

  return(results_df)
}
