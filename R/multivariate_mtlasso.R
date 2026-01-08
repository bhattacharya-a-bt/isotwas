#' Multivariate Multi-Task LASSO
#'
#' Fits a multi-task LASSO model for multivariate isoform prediction.
#' This method encourages joint feature selection across isoforms
#' using L21 regularization (row sparsity on the coefficient matrix).
#'
#' The L21 penalty enforces that the same SNPs are selected or excluded across
#' all isoforms, which is appropriate when SNPs are expected to have shared
#' effects on multiple isoforms of the same gene.
#'
#' @param X matrix, design matrix of SNP dosages (n x p)
#' @param Y matrix, matrix of isoform expression across columns (n x q)
#' @param regularization character, type of multi-task regularization:
#'   \itemize{
#'     \item "L21": L2,1 norm (row sparsity) - same SNPs selected across isoforms
#'     \item "Trace": trace norm (low-rank) - isoforms share latent structure (RMTL only)
#'     \item "Lasso": standard lasso applied jointly (RMTL only)
#'   }
#' @param lambda numeric or NULL, regularization parameter. If NULL, selected by CV.
#' @param lambda_seq numeric vector, sequence of lambda values to try in CV.
#'   If NULL, automatically generates a sequence.
#' @param nlambda int, number of lambda values if lambda_seq is NULL
#' @param lambda_min_ratio numeric, ratio of min to max lambda
#' @param Lam2 numeric, ridge penalty parameter (for RMTL backend). Default 0.
#' @param nfolds int, number of CV folds
#' @param standardize logical, standardize X before fitting. Default FALSE.
#' @param verbose logical, print progress
#' @param seed int, random seed
#' @param par logical, use parallel processing for CV folds
#' @param n.cores int, number of cores for parallel processing
#' @param backend character, implementation to use:
#'   \itemize{
#'     \item "fast": Custom fast implementation optimized for shared X (default)
#'     \item "rmtl": Use RMTL package (slower, but supports Trace/Lasso regularization)
#'   }
#'
#' @return isotwas_model object containing:
#'   \itemize{
#'     \item transcripts: list of transcript_model objects with weights, R2, pvalues
#'     \item best_lambda: optimal lambda from CV
#'     \item regularization: type of regularization used
#'   }
#'
#' @details
#' The L21 regularization solves:
#' \deqn{\min_W \frac{1}{2n}||Y - XW||_F^2 + \lambda ||W||_{2,1}}
#'
#' where \eqn{||W||_{2,1} = \sum_j ||W_j||_2} is the sum of L2 norms of rows,
#' encouraging entire rows (SNPs) to be zero or non-zero together.
#'
#' The fast backend exploits the shared X structure by precomputing X'X once,
#' making it much faster than general-purpose MTL solvers.
#'
#' @export
multivariate_mtlasso <- function(X,
                                  Y,
                                  regularization = c("L21", "Trace", "Lasso"),
                                  lambda = NULL,
                                  lambda_seq = NULL,
                                  nlambda = 20,
                                  lambda_min_ratio = 0.01,
                                  Lam2 = 0,
                                  nfolds = 5,
                                  standardize = FALSE,
                                  verbose = FALSE,
                                  seed = 123,
                                  par = FALSE,
                                  n.cores = 1,
                                  backend = c("fast", "rmtl")) {

  regularization <- match.arg(regularization)
  backend <- match.arg(backend)

  # Fast backend only supports L21

  if (backend == "fast" && regularization != "L21") {
    if (verbose) cat("Note: Fast backend only supports L21. Switching to RMTL.\n")
    backend <- "rmtl"
  }

  if (backend == "fast") {
    return(.mtlasso_fast(X, Y, lambda = lambda, lambda_seq = lambda_seq,
                         nlambda = nlambda, lambda_min_ratio = lambda_min_ratio,
                         nfolds = nfolds, standardize = standardize,
                         verbose = verbose, seed = seed,
                         par = par, n.cores = n.cores))
  } else {
    return(.mtlasso_rmtl(X, Y, regularization = regularization,
                         Lam1_seq = lambda_seq, Lam2 = Lam2,
                         nfolds = nfolds, verbose = verbose,
                         seed = seed, par = par, n.cores = n.cores))
  }
}


#' Fast Multi-Task LASSO with L21 penalty
#'
#' Custom implementation optimized for the case where all tasks share the same X.
#' Uses block coordinate descent with group soft-thresholding.
#'
#' @keywords internal
.mtlasso_fast <- function(X, Y, lambda = NULL, lambda_seq = NULL,
                          nlambda = 20, lambda_min_ratio = 0.01,
                          nfolds = 5, standardize = FALSE,
                          verbose = FALSE, seed = 123,
                          par = FALSE, n.cores = 1) {

  set.seed(seed)
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  tx_names <- colnames(Y)
  snp_names <- colnames(X)

  if (is.null(tx_names)) {
    tx_names <- paste0("Transcript", 1:q)
    colnames(Y) <- tx_names
  }
  if (is.null(snp_names)) {
    snp_names <- paste0("SNP", 1:p)
    colnames(X) <- snp_names
  }

  # Center Y
  Y_means <- colMeans(Y)
  Y_centered <- scale(Y, center = TRUE, scale = FALSE)

  # Standardize X if requested
  if (standardize) {
    X_means <- colMeans(X)
    X_sds <- apply(X, 2, sd)
    X_sds[X_sds == 0] <- 1
    X_scaled <- scale(X, center = TRUE, scale = TRUE)
  } else {
    X_means <- colMeans(X)
    X_scaled <- scale(X, center = TRUE, scale = FALSE)
    X_sds <- rep(1, p)
  }

  # Precompute X'X diagonal and X'Y (key optimization for shared X)
  XtX_diag <- colSums(X_scaled^2) / n
  XtY <- crossprod(X_scaled, Y_centered) / n

  # Compute lambda sequence
  if (is.null(lambda_seq)) {
    # Lambda_max: smallest lambda that gives all zeros
    # For L21, this is max over j of ||X_j'Y||_2 / n
    row_norms <- sqrt(rowSums(XtY^2))
    lambda_max <- max(row_norms)
    lambda_min <- lambda_max * lambda_min_ratio
    lambda_seq <- exp(seq(log(lambda_max), log(lambda_min), length.out = nlambda))
  }

  # If lambda provided directly, use it
  if (!is.null(lambda)) {
    best_lambda <- lambda
    if (verbose) cat(sprintf("Using provided lambda: %.6f\n", best_lambda))
  } else {
    # Cross-validation
    cv_folds <- caret::createFolds(1:n, k = nfolds, returnTrain = TRUE)

    if (verbose) cat("Running cross-validation...\n")

    # Define fold worker function
    run_fold <- function(fold_idx) {
      train_idx <- cv_folds[[fold_idx]]
      test_idx <- setdiff(1:n, train_idx)

      X_train <- X_scaled[train_idx, , drop = FALSE]
      Y_train <- Y_centered[train_idx, , drop = FALSE]
      X_test <- X_scaled[test_idx, , drop = FALSE]
      Y_test <- Y_centered[test_idx, , drop = FALSE]

      # Precompute for training data
      XtX_train <- colSums(X_train^2) / length(train_idx)
      XtY_train <- crossprod(X_train, Y_train) / length(train_idx)

      # Warm start path
      W <- matrix(0, p, q)
      fold_mse <- rep(0, length(lambda_seq))

      for (lam_idx in seq_along(lambda_seq)) {
        W <- .fit_mtl_l21(X_train, Y_train, lambda_seq[lam_idx],
                          XtX_diag = XtX_train, XtY = XtY_train,
                          W_init = W, max_iter = 200, tol = 1e-4)

        pred <- X_test %*% W
        fold_mse[lam_idx] <- mean((Y_test - pred)^2)
      }
      fold_mse
    }

    # Run CV folds (parallel or sequential)
    if (par && n.cores > 1) {
      fold_results <- parallel::mclapply(1:nfolds, run_fold, mc.cores = min(n.cores, nfolds))
    } else {
      fold_results <- lapply(1:nfolds, function(fold_idx) {
        if (verbose) cat(sprintf("  Fold %d/%d\n", fold_idx, nfolds))
        run_fold(fold_idx)
      })
    }

    # Aggregate MSE across folds
    cv_mse <- Reduce(`+`, fold_results) / nfolds
    best_idx <- which.min(cv_mse)
    best_lambda <- lambda_seq[best_idx]

    if (verbose) cat(sprintf("Best lambda: %.6f (index %d, CV MSE: %.4f)\n",
                             best_lambda, best_idx, cv_mse[best_idx]))
  }

  # Fit final model
  if (verbose) cat("Fitting final model...\n")
  W_final <- .fit_mtl_l21(X_scaled, Y_centered, best_lambda,
                          XtX_diag = XtX_diag, XtY = XtY,
                          W_init = NULL, max_iter = 1000, tol = 1e-5)

  # Rescale if standardized
  if (standardize) {
    W_final <- W_final / X_sds
  }

  # Get CV predictions for R² calculation (single pass)
  cv_folds <- caret::createFolds(1:n, k = nfolds, returnTrain = TRUE)
  cv_preds <- matrix(0, nrow = n, ncol = q)

  for (fold_idx in 1:nfolds) {
    train_idx <- cv_folds[[fold_idx]]
    test_idx <- setdiff(1:n, train_idx)

    X_train <- X_scaled[train_idx, , drop = FALSE]
    Y_train <- Y_centered[train_idx, , drop = FALSE]
    X_test <- X_scaled[test_idx, , drop = FALSE]

    XtX_train <- colSums(X_train^2) / length(train_idx)
    XtY_train <- crossprod(X_train, Y_train) / length(train_idx)

    W_fold <- .fit_mtl_l21(X_train, Y_train, best_lambda,
                           XtX_diag = XtX_train, XtY = XtY_train,
                           W_init = NULL, max_iter = 500, tol = 1e-4)
    cv_preds[test_idx, ] <- X_test %*% W_fold
  }

  # Add back means
  final_preds <- cv_preds + matrix(Y_means, nrow = n, ncol = q, byrow = TRUE)

  # Compute R² and p-values
  r2_vec <- numeric(q)
  p_vec <- numeric(q)

  for (j in 1:q) {
    y <- Y[, j]
    pred <- final_preds[, j]

    ss_res <- sum((y - pred)^2)
    ss_tot <- sum((y - mean(y))^2)

    if (ss_tot > 0) {
      r2 <- 1 - ss_res / ss_tot
      r2_adj <- 1 - (1 - r2) * (n - 1) / (n - 2)
      r2_vec[j] <- max(0, r2_adj)
    }

    if (sd(pred) > 0) {
      r <- cor(y, pred)
      t_stat <- r * sqrt((n - 2) / (1 - r^2))
      p_vec[j] <- 2 * pt(-abs(t_stat), df = n - 2)
    } else {
      p_vec[j] <- 1
    }
  }

  # Build output
  transcripts <- list()
  for (j in 1:q) {
    tx_name <- tx_names[j]
    weights_df <- tibble::tibble(
      SNP = snp_names,
      Weight = W_final[, j]
    )
    weights_df <- subset(weights_df, Weight != 0)

    transcripts[[tx_name]] <- create_transcript_model(
      transcript_id = tx_name,
      weights = weights_df,
      r2 = r2_vec[j],
      pvalue = p_vec[j],
      predicted = final_preds[, j]
    )
  }

  result <- create_isotwas_model(
    method = "mtlasso_l21",
    transcripts = transcripts,
    n_samples = n,
    n_snps = p
  )

  result$regularization <- "L21"
  result$best_lambda <- best_lambda
  result$lambda_seq <- lambda_seq
  result$backend <- "fast"
  if (exists("cv_mse")) result$cv_mse <- cv_mse

  return(result)
}


#' Fit L21-penalized multi-task regression via block coordinate descent
#'
#' @keywords internal
.fit_mtl_l21 <- function(X, Y, lambda, XtX_diag = NULL, XtY = NULL,
                         W_init = NULL, max_iter = 500, tol = 1e-4) {
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)

  # Initialize
  if (is.null(W_init)) {
    W <- matrix(0, p, q)
  } else {
    W <- W_init
  }

  # Precompute if not provided
  if (is.null(XtX_diag)) {
    XtX_diag <- colSums(X^2) / n
  }
  if (is.null(XtY)) {
    XtY <- crossprod(X, Y) / n
  }

  # Residual
  R <- Y - X %*% W

  for (iter in 1:max_iter) {
    W_old <- W
    max_change <- 0

    # Block coordinate descent over SNPs (rows of W)
    for (j in 1:p) {
      w_j <- W[j, ]

      # Add back contribution
      R <- R + outer(X[, j], w_j)

      # Gradient: (1/n) X_j' R
      grad_j <- as.vector(crossprod(X[, j], R)) / n

      # Group soft-thresholding for L21
      # Update: w_j = S_λ(grad_j / XtX_jj) where S_λ is group soft-threshold
      # S_λ(v) = v * max(0, 1 - λ/||v||_2)

      grad_norm <- sqrt(sum(grad_j^2))

      if (grad_norm > lambda) {
        # Non-zero update
        shrink <- 1 - lambda / grad_norm
        w_j_new <- shrink * grad_j / XtX_diag[j]
      } else {
        # Zero out entire row
        w_j_new <- rep(0, q)
      }

      W[j, ] <- w_j_new
      R <- R - outer(X[, j], w_j_new)

      max_change <- max(max_change, max(abs(w_j_new - w_j)))
    }

    if (max_change < tol) {
      break
    }
  }

  return(W)
}


#' RMTL backend for Multi-Task LASSO
#'
#' Uses the RMTL package. Slower but supports Trace and Lasso regularization.
#'
#' @keywords internal
.mtlasso_rmtl <- function(X, Y, regularization = "L21",
                          Lam1_seq = NULL, Lam2 = 0,
                          nfolds = 5, verbose = FALSE, seed = 123,
                          par = FALSE, n.cores = 1) {

  set.seed(seed)
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  tx_names <- colnames(Y)
  snp_names <- colnames(X)

  if (is.null(tx_names)) {
    tx_names <- paste0("Transcript", 1:q)
    colnames(Y) <- tx_names
  }
  if (is.null(snp_names)) {
    snp_names <- paste0("SNP", 1:p)
    colnames(X) <- snp_names
  }

  # RMTL format: list of X matrices and Y vectors
  X_list <- lapply(1:q, function(i) as.matrix(X))
  Y_list <- lapply(1:q, function(i) Y[, i])

  if (is.null(Lam1_seq)) {
    Lam1_seq <- 10^seq(1, -3, length.out = 15)
  }

  if (verbose) cat("Running RMTL cross-validation...\n")

  cv_result <- tryCatch({
    RMTL::cvMTL(
      X = X_list, Y = Y_list,
      type = "Regression",
      Regularization = regularization,
      Lam1_seq = Lam1_seq,
      Lam2 = Lam2,
      opts = list(init = 0, tol = 1e-4, maxIter = 500),
      nfolds = nfolds,
      stratify = FALSE,
      parallel = par,
      ncores = n.cores
    )
  }, error = function(e) {
    stop("RMTL failed: ", e$message)
  })

  best_lambda <- cv_result$Lam1.min
  if (verbose) cat(sprintf("Best lambda: %.6f\n", best_lambda))

  # Fit final model
  if (verbose) cat("Fitting final model...\n")
  final_model <- RMTL::MTL(
    X = X_list, Y = Y_list,
    type = "Regression",
    Regularization = regularization,
    Lam1 = best_lambda,
    Lam2 = Lam2,
    opts = list(init = 0, tol = 1e-5, maxIter = 1000)
  )

  W <- final_model$W
  intercepts <- final_model$C

  # CV predictions using RMTL predict (faster than refitting)
  cv_folds <- caret::createFolds(1:n, k = nfolds, returnTrain = TRUE)
  cv_preds <- matrix(0, nrow = n, ncol = q)

  if (verbose) cat("Computing CV predictions...\n")

  for (fold_idx in 1:nfolds) {
    train_idx <- cv_folds[[fold_idx]]
    test_idx <- setdiff(1:n, train_idx)

    X_train_list <- lapply(1:q, function(i) as.matrix(X[train_idx, , drop = FALSE]))
    Y_train_list <- lapply(1:q, function(i) Y[train_idx, i])
    X_test <- X[test_idx, , drop = FALSE]

    fold_model <- tryCatch({
      RMTL::MTL(
        X = X_train_list, Y = Y_train_list,
        type = "Regression",
        Regularization = regularization,
        Lam1 = best_lambda,
        Lam2 = Lam2,
        opts = list(init = 0, tol = 1e-4, maxIter = 500)
      )
    }, error = function(e) NULL)

    if (!is.null(fold_model)) {
      cv_preds[test_idx, ] <- as.matrix(X_test) %*% fold_model$W +
        matrix(fold_model$C, nrow = length(test_idx), ncol = q, byrow = TRUE)
    }
  }

  # R² and p-values
  r2_vec <- numeric(q)
  p_vec <- numeric(q)

  for (j in 1:q) {
    y <- Y[, j]
    pred <- cv_preds[, j]

    ss_res <- sum((y - pred)^2)
    ss_tot <- sum((y - mean(y))^2)

    if (ss_tot > 0) {
      r2 <- 1 - ss_res / ss_tot
      r2_adj <- 1 - (1 - r2) * (n - 1) / (n - 2)
      r2_vec[j] <- max(0, r2_adj)
    }

    if (sd(pred) > 0) {
      r <- cor(y, pred)
      t_stat <- r * sqrt((n - 2) / (1 - r^2))
      p_vec[j] <- 2 * pt(-abs(t_stat), df = n - 2)
    } else {
      p_vec[j] <- 1
    }
  }

  # Build output
  transcripts <- list()
  for (j in 1:q) {
    tx_name <- tx_names[j]
    weights_df <- tibble::tibble(
      SNP = snp_names,
      Weight = W[, j]
    )
    weights_df <- subset(weights_df, Weight != 0)

    transcripts[[tx_name]] <- create_transcript_model(
      transcript_id = tx_name,
      weights = weights_df,
      r2 = r2_vec[j],
      pvalue = p_vec[j],
      predicted = cv_preds[, j]
    )
  }

  result <- create_isotwas_model(
    method = paste0("mtlasso_", tolower(regularization)),
    transcripts = transcripts,
    n_samples = n,
    n_snps = p
  )

  result$regularization <- regularization
  result$best_lambda <- best_lambda
  result$Lam2 <- Lam2
  result$intercepts <- intercepts
  result$backend <- "rmtl"

  return(result)
}


#' Multivariate Multi-Task LASSO with multiple regularization comparison
#'
#' Runs MTLASSO with different regularization types and returns the best model
#' based on cross-validated MSE.
#'
#' @param X matrix, design matrix of SNP dosages (n x p)
#' @param Y matrix, matrix of isoform expression across columns (n x q)
#' @param regularizations character vector, regularization types to compare
#' @param nfolds int, number of CV folds
#' @param verbose logical, print progress
#' @param seed int, random seed
#' @param par logical, use parallel processing for CV folds
#' @param n.cores int, number of cores for parallel processing
#'
#' @return isotwas_model object with best regularization
#'
#' @export
multivariate_mtlasso_cv <- function(X,
                                     Y,
                                     regularizations = c("L21", "Trace", "Lasso"),
                                     nfolds = 5,
                                     verbose = FALSE,
                                     seed = 123,
                                     par = FALSE,
                                     n.cores = 1) {

  best_mse <- Inf
  best_model <- NULL

  for (reg in regularizations) {
    if (verbose) cat(sprintf("Trying regularization: %s\n", reg))

    # Use fast backend for L21, RMTL for others
    backend <- if (reg == "L21") "fast" else "rmtl"

    model <- tryCatch({
      multivariate_mtlasso(
        X = X, Y = Y,
        regularization = reg,
        nfolds = nfolds,
        verbose = FALSE,
        seed = seed,
        par = par,
        n.cores = n.cores,
        backend = backend
      )
    }, error = function(e) {
      if (verbose) cat(sprintf("  Failed: %s\n", e$message))
      NULL
    })

    if (!is.null(model)) {
      mse <- mean(sapply(model$transcripts, function(tx) {
        y <- Y[, tx$transcript_id]
        pred <- tx$predicted
        mean((y - pred)^2)
      }))

      if (verbose) cat(sprintf("  Mean MSE: %.4f, Mean R2: %.4f\n",
                               mse, mean(model$summary$r2)))

      if (mse < best_mse) {
        best_mse <- mse
        best_model <- model
      }
    }
  }

  if (is.null(best_model)) {
    stop("All regularization types failed")
  }

  if (verbose) cat(sprintf("Best regularization: %s\n", best_model$regularization))

  return(best_model)
}
