#' isoTWAS Model Classes and Methods
#'
#' S3 classes for representing isoTWAS prediction model outputs with
#' intuitive structure and accessor methods.

#' Create an isotwas_model object
#'
#' @param method character, name of the prediction method
#' @param transcripts named list of transcript model results
#' @param n_samples int, number of samples used for training
#' @param n_snps int, total number of SNPs in input
#'
#' @return An isotwas_model object
#'
#' @export
create_isotwas_model <- function(method,
                                  transcripts,
                                  n_samples,
                                  n_snps) {

  # Ensure transcripts is a named list
  if (is.null(names(transcripts))) {
    stop("transcripts must be a named list")
  }

  # Build summary data frame
  summary_df <- data.frame(
    transcript = names(transcripts),
    r2 = sapply(transcripts, function(x) x$r2),
    pvalue = sapply(transcripts, function(x) x$pvalue),
    n_snps_selected = sapply(transcripts, function(x) nrow(x$weights)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  structure(
    list(
      method = method,
      transcripts = transcripts,
      summary = summary_df,
      n_transcripts = length(transcripts),
      n_samples = n_samples,
      n_snps_total = n_snps
    ),
    class = "isotwas_model"
  )
}

#' Create a single transcript model entry
#'
#' Helper function to create a properly structured transcript model
#'
#' @param transcript_id character, transcript/isoform identifier
#' @param weights tibble/data.frame with SNP and Weight columns
#' @param r2 numeric, cross-validated R-squared
#' @param pvalue numeric, p-value from prediction correlation
#' @param predicted numeric vector, cross-validated predictions
#'
#' @return A list representing a single transcript model
#'
#' @export
create_transcript_model <- function(transcript_id,
                                     weights,
                                     r2,
                                     pvalue,
                                     predicted) {
  list(
    transcript_id = transcript_id,
    weights = weights,
    r2 = r2,
    pvalue = pvalue,
    predicted = predicted
  )
}

#' Print method for isotwas_model
#'
#' @param x An isotwas_model object
#' @param ... Additional arguments (unused)
#'
#' @export
print.isotwas_model <- function(x, ...) {
  cat("isoTWAS Prediction Model\n")
  cat("========================\n")
  cat(sprintf("Method: %s\n", x$method))
  cat(sprintf("Transcripts: %d\n", x$n_transcripts))
  cat(sprintf("Samples: %d\n", x$n_samples))
  cat(sprintf("Total SNPs: %d\n", x$n_snps_total))
  cat("\n")
  cat("Performance Summary:\n")
  cat(sprintf("  Mean R2: %.4f\n", mean(x$summary$r2, na.rm = TRUE)))
  cat(sprintf("  Median R2: %.4f\n", median(x$summary$r2, na.rm = TRUE)))
  cat(sprintf("  R2 range: [%.4f, %.4f]\n",
              min(x$summary$r2, na.rm = TRUE),
              max(x$summary$r2, na.rm = TRUE)))
  cat(sprintf("  Mean SNPs selected: %.1f\n",
              mean(x$summary$n_snps_selected, na.rm = TRUE)))
  cat("\n")
  cat("Use summary() for per-transcript details\n")
  cat("Access transcripts with $transcripts$<name> or [[\"<name>\"]]\n")
  invisible(x)
}

#' Summary method for isotwas_model
#'
#' @param object An isotwas_model object
#' @param ... Additional arguments (unused)
#'
#' @return A data.frame with per-transcript summary statistics
#'
#' @export
summary.isotwas_model <- function(object, ...) {
  cat(sprintf("isoTWAS Model Summary: %s\n", object$method))
  cat(paste(rep("-", 50), collapse = ""), "\n")
  print(object$summary, row.names = FALSE)
  invisible(object$summary)
}

#' Extract a transcript model from isotwas_model
#'
#' @param x An isotwas_model object
#' @param name Transcript name to extract
#'
#' @export
`$.isotwas_model` <- function(x, name) {
  if (name %in% names(x)) {
    x[[name]]
  } else if (name %in% names(x$transcripts)) {
    x$transcripts[[name]]
  } else {
    NULL
  }
}

#' Subset isotwas_model by transcript name
#'
#' @param x An isotwas_model object
#' @param i Transcript name(s) or index
#'
#' @export
`[.isotwas_model` <- function(x, i) {
  if (is.character(i)) {
    x$transcripts[i]
  } else {
    x$transcripts[i]
  }
}

#' Get weights for all transcripts as a matrix
#'
#' @param model An isotwas_model object
#' @param snps Optional character vector of SNP names to include (default: all)
#'
#' @return A matrix with SNPs as rows and transcripts as columns
#'
#' @export
get_weight_matrix <- function(model, snps = NULL) {
  if (!inherits(model, "isotwas_model")) {
    stop("model must be an isotwas_model object")
  }

  # Get all unique SNPs across transcripts
  all_snps <- unique(unlist(lapply(model$transcripts, function(x) x$weights$SNP)))

  if (!is.null(snps)) {
    all_snps <- intersect(all_snps, snps)
  }

  # Build weight matrix
  weight_mat <- matrix(0,
                       nrow = length(all_snps),
                       ncol = model$n_transcripts,
                       dimnames = list(all_snps, names(model$transcripts)))

  for (tx_name in names(model$transcripts)) {
    tx <- model$transcripts[[tx_name]]
    matched <- match(tx$weights$SNP, all_snps)
    valid <- !is.na(matched)
    weight_mat[matched[valid], tx_name] <- tx$weights$Weight[valid]
  }

  weight_mat
}

#' Get predicted values as a matrix
#'
#' @param model An isotwas_model object
#'
#' @return A matrix with samples as rows and transcripts as columns
#'
#' @export
get_prediction_matrix <- function(model) {
  if (!inherits(model, "isotwas_model")) {
    stop("model must be an isotwas_model object")
  }

  pred_mat <- sapply(model$transcripts, function(x) x$predicted)
  colnames(pred_mat) <- names(model$transcripts)
  pred_mat
}


# === isotwas_result class for compute_isotwas output ===

#' Create an isotwas_result object
#'
#' @param best_models Named list of best models per transcript
#' @param method_results List of isotwas_model objects (one per method)
#' @param comparison_df Data frame comparing methods
#' @param tx2gene Transcript to gene coefficients
#'
#' @return An isotwas_result object
#'
#' @export
create_isotwas_result <- function(best_models,
                                   method_results = NULL,
                                   comparison_df = NULL,
                                   tx2gene = NULL) {

  structure(
    list(
      best_models = best_models,
      method_results = method_results,
      comparison = comparison_df,
      tx2gene = tx2gene,
      n_transcripts = length(best_models$transcripts)
    ),
    class = "isotwas_result"
  )
}

#' Print method for isotwas_result
#'
#' @param x An isotwas_result object
#' @param ... Additional arguments (unused)
#'
#' @export
print.isotwas_result <- function(x, ...) {
  cat("isoTWAS Analysis Result\n")
  cat("=======================\n")
  cat(sprintf("Transcripts modeled: %d\n", x$n_transcripts))

  if (!is.null(x$method_results)) {
    cat(sprintf("Methods compared: %s\n",
                paste(names(x$method_results), collapse = ", ")))
  }

  cat("\nBest Model Performance:\n")
  print(x$best_models)

  if (!is.null(x$comparison) && nrow(x$comparison) > 0) {
    cat("\nMethod Comparison (Mean R2):\n")
    method_means <- colMeans(x$comparison[, -1, drop = FALSE], na.rm = TRUE)
    for (m in names(method_means)) {
      cat(sprintf("  %s: %.4f\n", m, method_means[m]))
    }
  }

  if (!is.null(x$tx2gene) && is.data.frame(x$tx2gene)) {
    cat(sprintf("\nTranscript-to-gene model R2: %.4f\n", x$tx2gene$R2[1]))
  }

  invisible(x)
}

#' Summary method for isotwas_result
#'
#' @param object An isotwas_result object
#' @param ... Additional arguments (unused)
#'
#' @export
summary.isotwas_result <- function(object, ...) {
  cat("isoTWAS Result Summary\n")
  cat(paste(rep("=", 50), collapse = ""), "\n\n")

  cat("Best Models per Transcript:\n")
  summary(object$best_models)

  if (!is.null(object$comparison)) {
    cat("\n\nR2 Comparison Across Methods:\n")
    print(object$comparison, row.names = FALSE)
  }

  invisible(object)
}

#' Convert legacy model list to isotwas_model
#'
#' Converts the old list-of-lists format to the new isotwas_model class
#'
#' @param model_list List in legacy format (list of lists with Transcript, Model, R2, P, Pred)
#' @param method Character string identifying the method
#' @param n_samples Number of samples
#' @param n_snps Total number of SNPs
#'
#' @return An isotwas_model object
#'
#' @export
convert_legacy_model <- function(model_list, method, n_samples, n_snps) {
  transcripts <- list()

  for (item in model_list) {
    tx_name <- item$Transcript
    transcripts[[tx_name]] <- create_transcript_model(
      transcript_id = tx_name,
      weights = item$Model,
      r2 = item$R2,
      pvalue = item$P,
      predicted = item$Pred
    )
  }

  create_isotwas_model(
    method = method,
    transcripts = transcripts,
    n_samples = n_samples,
    n_snps = n_snps
  )
}
