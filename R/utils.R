"%^%" <- function(x, n) {
    with(eigen(x), vectors %*% (values^n * t(vectors)))
    }

.onUnload <- function (libpath) { library.dynam.unload("isotwas", libpath)}

calc.mse <- function(obs, pred){
    if(is.vector(obs)) obs <- as.matrix(obs)
    if(is.vector(pred)) pred <- as.matrix(pred)

    n <- nrow(obs)
    rss <- colSums((obs - pred)^2, na.rm = TRUE)
    rss/n
}

calc.r2 <- function(i,obs,pred){

    r2 = summary(stats::lm(obs[,i]~pred[,i]))$adj.r.sq
    ct = stats::cor.test(obs[,i],pred[,i])
    if (is.na(r2)){r2 = 0}
    p = ct$p.value
    return(list(R2 = r2,
                P = p))

}


get_best <- function(list_mods, G = G){
    # Handle new isotwas_model objects
    if (length(list_mods) > 0 && inherits(list_mods[[1]], "isotwas_model")) {
        return(get_best_isotwas(list_mods, G))
    }

    # Legacy behavior for old list format
    r2 = matrix(unlist(sapply(list_mods, function(y) sapply(y,function(x) x$R2))),
                nrow = G)
    l_nrows = sapply(list_mods, function(y) sapply(y,function(x) nrow(x$Model)))
    r2[l_nrows == 0] = 0
    r2[is.na(r2)] = 0
    bundle = cbind(apply(r2,1,which.max),1:G)
    return(lapply(1:G,
                  function(i) list_mods[[bundle[i,1]]][[bundle[i,2]]]))
}

#' Select best model per transcript from multiple isotwas_model objects
#'
#' @param model_list List of isotwas_model objects
#' @param G Number of transcripts
#'
#' @return An isotwas_model with the best model per transcript
get_best_isotwas <- function(model_list, G) {
    # Get transcript names from first model
    tx_names <- names(model_list[[1]]$transcripts)

    # Build R2 matrix: rows = transcripts, cols = methods
    r2_mat <- matrix(nrow = G, ncol = length(model_list))
    method_names <- sapply(model_list, function(m) m$method)
    colnames(r2_mat) <- method_names
    rownames(r2_mat) <- tx_names

    for (j in seq_along(model_list)) {
        mod <- model_list[[j]]
        for (i in seq_along(tx_names)) {
            tx <- tx_names[i]
            if (tx %in% names(mod$transcripts)) {
                r2_val <- mod$transcripts[[tx]]$r2
                n_snps <- nrow(mod$transcripts[[tx]]$weights)
                # Set R2 to 0 if no SNPs selected
                if (is.null(n_snps) || n_snps == 0) r2_val <- 0
                r2_mat[i, j] <- ifelse(is.na(r2_val), 0, r2_val)
            } else {
                r2_mat[i, j] <- 0
            }
        }
    }

    # Select best method per transcript
    best_idx <- apply(r2_mat, 1, which.max)

    # Build combined model with best per transcript
    best_transcripts <- list()
    best_methods <- character(G)

    for (i in seq_along(tx_names)) {
        tx <- tx_names[i]
        best_j <- best_idx[i]
        best_transcripts[[tx]] <- model_list[[best_j]]$transcripts[[tx]]
        best_transcripts[[tx]]$method <- method_names[best_j]
        best_methods[i] <- method_names[best_j]
    }

    # Get sample/SNP counts from first model
    n_samples <- model_list[[1]]$n_samples
    n_snps <- model_list[[1]]$n_snps_total

    result <- create_isotwas_model(
        method = "best",
        transcripts = best_transcripts,
        n_samples = n_samples,
        n_snps = n_snps
    )

    # Add comparison info
    result$r2_comparison <- r2_mat
    result$best_methods <- data.frame(
        transcript = tx_names,
        best_method = best_methods,
        r2 = sapply(seq_along(tx_names), function(i) r2_mat[i, best_idx[i]]),
        stringsAsFactors = FALSE
    )

    return(result)
}

cluster_weight = function(x){
  if (length(x) == 0){return(list(sum = 0,
                                  index = NA))}
  if (length(x) == 1){return(list(sum = 0,
                                  index = 1))}
  p = stats::p.adjust(stats::pnorm(
    as.numeric(scale(x)),
    lower.tail=F),
    'fdr')
  if (sum(p < 0.05) > 0){
    return(list(sum = sum(p < 0.05),
                index = which(p < 0.05)))
  }
  if (sum(p < .25) > 0){
    return(list(sum = sum(p < 0.25),
                index = which(p < 0.25)))
  }
  if (sum(p < .25) == 0){
    return(list(sum = length(x),
                index = 1:length(x)))
  }
}

