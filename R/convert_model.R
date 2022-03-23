#' Convert isotwas model list to tibble
#'
#' The function converts an isotwas model and returns a tibble
#'
#' @param model, isotwas model
#' @param snp_annot, annotations for SNPs, typically map obj from bigsnpr
#' @param snp_var, column name of SNP marker id
#'
#' @return tibble of the model
#'
#' @export
convert_model = function(model,
                         snp_annot = NULL,
                         snp_var = NULL){

  mod = as.data.frame(matrix(nrow = 0,
                             ncol = 5))
  colnames(mod) = c('SNP','Weight','Transcript','R2','P')

  for (i in 1:length(model$Model)){

    this.model = as.data.frame(model$Model[[i]]$Model)
    this.model$Transcript = model$Model[[i]]$Transcript
    this.model$R2 = model$Model[[i]]$R2
    this.model$R2.P = model$Model[[i]]$P
    mod = rbind(mod,
                this.model)

  }

  if (!is.null(snp_annot)){
    if (!'SNP' %in% colnames(snp_annot)){
      colnames(snp_annot)[which(colnames(snp_annot) == snp_var)] = 'SNP'
    }
    mod = merge(mod,snp_annot,by='SNP')

  }

  return(tibble::as_tibble(mod))

}
