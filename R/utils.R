"%^%" <- function(x, n) {
    with(eigen(x), vectors %*% (values^n * t(vectors)))
    }

.onUnload <- function (libpath) { library.dynam.unload("isotwas", libpath)}

#' Methods compute_isotwas can run, in comparison-column order
#'
#' Each fitting function labels its output $method with its name here.
#'
#' @return character vector of method names
#'
#' @noRd
isotwas_methods <- function() {
    c('mrce_lasso', 'multi_enet', 'joinet', 'spls',
      'sgl', 'mtlasso', 'stacking', 'graph_reg', 'univariate')
}

#' mclapply with reproducible, core-count-independent per-job seeds
#'
#' Each job gets its own seed drawn from \code{seed}, so results do not depend
#' on mc.cores (RNGkind("L'Ecuyer-CMRG") alone would). The caller's RNG state
#' is restored on exit.
#'
#' @param X vector or list to iterate over
#' @param FUN function applied to each element
#' @param ... further arguments to FUN
#' @param mc.cores int, passed to parallel::mclapply
#' @param seed int, the seed the per-job seeds are derived from
#'
#' @return list, as parallel::mclapply
#'
#' @noRd
.mclapply_seeded <- function(X, FUN, ..., mc.cores = 1L, seed) {

    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
        old_seed <- get(".Random.seed", envir = globalenv())
        on.exit(assign(".Random.seed", old_seed, envir = globalenv()),
                add = TRUE)
    }

    set.seed(seed)
    job_seeds <- sample.int(.Machine$integer.max, length(X))

    parallel::mclapply(seq_along(X), function(k) {
        set.seed(job_seeds[k])
        FUN(X[[k]], ...)
    }, mc.cores = mc.cores)
}

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


#' Score a prediction against an observation
#'
#' "refit" is calc.r2(): the adjusted R-squared of lm(obs ~ pred). "oos" is
#' the adjusted out-of-sample R-squared 1 - ss_res/ss_tot, which mrce_lasso,
#' sgl, mtlasso, stacking and graph_reg compute; it is negative when a
#' prediction is worse than the mean.
#'
#' @param obs numeric vector, observed values
#' @param pred numeric vector, predicted values
#' @param method character, "refit" or "oos"
#'
#' @return list(R2, P)
#'
#' @noRd
.score_r2 <- function(obs, pred, method = c("refit", "oos")) {

    method <- match.arg(method)
    n <- length(obs)

    if (method == "refit") {
        r2 <- summary(stats::lm(obs ~ pred))$adj.r.sq
        if (is.na(r2)) r2 <- 0
    } else {
        ss_res <- sum((obs - pred)^2)
        ss_tot <- sum((obs - mean(obs))^2)
        if (ss_tot > 0) {
            r2_raw <- 1 - ss_res / ss_tot
            r2 <- 1 - (1 - r2_raw) * (n - 1) / (n - 2)
        } else {
            r2 <- 0
        }
    }

    p <- if (stats::sd(pred) > 0 && stats::sd(obs) > 0) {
        stats::cor.test(obs, pred)$p.value
    } else {
        1
    }

    list(R2 = r2, P = p)
}


get_best <- function(list_mods, G = G, empty_r2 = 0){
    # Handle new isotwas_model objects
    if (length(list_mods) > 0 && inherits(list_mods[[1]], "isotwas_model")) {
        return(get_best_isotwas(list_mods, G, empty_r2 = empty_r2))
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

#' Resolve a transcript's method label through a nested ensemble
#'
#' A nested ensemble (the univariate entry) carries a $best_methods table;
#' its leaf label is used in place of the parent's "best".
#'
#' @param mod an isotwas_model, possibly carrying a $best_methods table
#' @param tx character, transcript id
#' @param fallback character, the label to use when mod is a flat method
#'
#' @return character, the leaf method name
#'
#' @noRd
.leaf_method <- function(mod, tx, fallback) {
    bm <- mod$best_methods
    if (is.null(bm) || !is.data.frame(bm) || !"transcript" %in% names(bm)) {
        return(fallback)
    }
    hit <- match(tx, bm$transcript)
    if (is.na(hit)) return(fallback)
    leaf <- bm$best_method[hit]
    if (is.na(leaf) || !nzchar(leaf)) fallback else leaf
}

#' Select best model per transcript from multiple isotwas_model objects
#'
#' @param model_list List of isotwas_model objects
#' @param G Number of transcripts
#' @param empty_r2 numeric, the score given to a model that selected no SNPs,
#'   or to a transcript absent from a model. 0 is the historical value; with
#'   -Inf, a transcript for which every model is empty is labelled NA.
#'
#' @return An isotwas_model with the best model per transcript
get_best_isotwas <- function(model_list, G, empty_r2 = 0) {
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
                if (is.null(n_snps) || n_snps == 0) r2_val <- empty_r2
                r2_mat[i, j] <- ifelse(is.na(r2_val), empty_r2, r2_val)
            } else {
                r2_mat[i, j] <- empty_r2
            }
        }
    }

    # Select best method per transcript. which.max() returns 1 for a row with
    # no finite score, so those transcripts are flagged and labelled NA.
    unfittable <- apply(r2_mat, 1, function(x) all(!is.finite(x)))
    best_idx <- apply(r2_mat, 1, which.max)

    # Build combined model with best per transcript
    best_transcripts <- list()
    best_methods <- character(G)

    for (i in seq_along(tx_names)) {
        tx <- tx_names[i]
        best_j <- best_idx[i]
        best_transcripts[[tx]] <- model_list[[best_j]]$transcripts[[tx]]
        leaf <- .leaf_method(model_list[[best_j]], tx, method_names[best_j])
        if (unfittable[i]) leaf <- NA_character_
        best_transcripts[[tx]]$method <- leaf
        best_methods[i] <- leaf
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

