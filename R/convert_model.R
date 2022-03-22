#' Convert isotwas model list to tibble
#'
#' The function converts an isotwas model from compute_isotwas() and returns a tibble
#'
#' @param model
#'
#' @return tibble of the model
#'
#' @export
convert_model = function(model){

  mod = as.data.frame(matrix(nrow = 0,
                             ncol = 5))
  colnames(mod) = c('SNP','Weight','Transcript','R2','P')

  for (i in 1:length(model)){

    this.model = as.data.frame(model$Model[[i]]$Model)
    this.model$Transcript = model$Model[[i]]$Transcript
    this.model$R2 = model$Model[[i]]$R2
    this.model$R2.P = model$Model[[i]]$P
    mod = rbind(mod,
                this.model)

  }

  return(tibble::as_tibble(mod))

}
