#' Multivariate Sparse Group LASSO
#'
#' Fits a sparse group LASSO model for multivariate isoform prediction.
#' This method encourages both:
#' - Group-level sparsity: same SNPs tend to affect multiple isoforms (shared effects)
#' - Within-group sparsity: isoform-specific effects
#'
#' The penalty is: alpha * ||B||_1 + (1-alpha) * sum_j ||B_j||_2
#' where B_j is the j-th row of B (effects of SNP j across all isoforms)
#'
#' @param X matrix, design matrix of SNP dosages
#' @param Y matrix, matrix of G isoform expression across columns
#' @param alpha numeric, mixing parameter between group and individual sparsity (0-1)
#'   alpha=1 is pure LASSO, alpha=0 is pure group LASSO
#' @param nlambda int, number of lambda values to try
#' @param lambda_min_ratio numeric, ratio of lambda_min to lambda_max
#' @param nfolds int, number of CV folds
#' @param standardize logical, standardize X before fitting
#' @param verbose logical
#' @param seed int, random seed
#' @param par logical, use parallel processing for CV folds
#' @param n.cores int, number of cores for parallel processing
#'
#' @return isotwas_model object
#'
#' @export
multivariate_sgl <- function(X,
                              Y,
                              alpha = 0.5,
                              nlambda = 20,
                              lambda_min_ratio = 0.01,
                              nfolds = 5,
                              standardize = FALSE,
                              verbose = FALSE,
                              seed = 123,
                              par = FALSE,
                              n.cores = 1) {

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

  # Optionally standardize X
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

  # Compute lambda sequence
  # Lambda_max: smallest lambda that gives all zeros
  lambda_max <- .compute_lambda_max_sgl(X_scaled, Y_centered, alpha)
  lambda_min <- lambda_max * lambda_min_ratio
  lambda_seq <- exp(seq(log(lambda_max), log(lambda_min), length.out = nlambda))

  # Create CV folds
  cv_folds <- caret::createFolds(1:n, k = nfolds, returnTrain = TRUE)

  # CV to find best lambda
  if (verbose) cat("Running cross-validation for lambda selection...\n")

  # Define fold worker function
  run_fold <- function(fold_idx) {
    train_idx <- cv_folds[[fold_idx]]
    test_idx <- setdiff(1:n, train_idx)

    X_train <- X_scaled[train_idx, , drop = FALSE]
    Y_train <- Y_centered[train_idx, , drop = FALSE]
    X_test <- X_scaled[test_idx, , drop = FALSE]
    Y_test <- Y_centered[test_idx, , drop = FALSE]

    # Warm start: use previous lambda solution
    B_warm <- matrix(0, nrow = p, ncol = q)
    fold_mse <- matrix(0, nrow = nlambda, ncol = q)

    for (lam_idx in 1:nlambda) {
      B_warm <- .fit_sgl(X_train, Y_train, lambda_seq[lam_idx], alpha,
                         B_init = B_warm, max_iter = 500, tol = 1e-4)

      # Predictions
      pred <- X_test %*% B_warm
      fold_mse[lam_idx, ] <- colMeans((Y_test - pred)^2)
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

  # Select best lambda (minimize mean MSE across transcripts)
  mean_mse <- rowMeans(cv_mse)
  best_lambda_idx <- which.min(mean_mse)
  best_lambda <- lambda_seq[best_lambda_idx]

  if (verbose) cat(sprintf("Best lambda: %.6f (index %d)\n", best_lambda, best_lambda_idx))

  # Fit final model on full data
  if (verbose) cat("Fitting final model...\n")
  B_final <- .fit_sgl(X_scaled, Y_centered, best_lambda, alpha,
                      B_init = NULL, max_iter = 1000, tol = 1e-5)

  # Rescale coefficients if standardized
  if (standardize) {
    B_final <- B_final / X_sds
  }

  # Get CV predictions at best lambda for R² calculation
  cv_preds <- matrix(0, nrow = n, ncol = q)
  for (fold_idx in 1:nfolds) {
    train_idx <- cv_folds[[fold_idx]]
    test_idx <- setdiff(1:n, train_idx)

    X_train <- X_scaled[train_idx, , drop = FALSE]
    Y_train <- Y_centered[train_idx, , drop = FALSE]
    X_test <- X_scaled[test_idx, , drop = FALSE]

    B_fold <- .fit_sgl(X_train, Y_train, best_lambda, alpha,
                       B_init = NULL, max_iter = 500, tol = 1e-4)
    cv_preds[test_idx, ] <- X_test %*% B_fold
  }

  # Add back means for final predictions
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
      Weight = B_final[, j]
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
    method = "sparse_group_lasso",
    transcripts = transcripts,
    n_samples = n,
    n_snps = p
  )

  # Add SGL-specific info
  result$alpha <- alpha
  result$best_lambda <- best_lambda
  result$cv_mse <- cv_mse
  result$lambda_seq <- lambda_seq

  return(result)
}


#' Compute lambda_max for Sparse Group LASSO
#' @keywords internal
.compute_lambda_max_sgl <- function(X, Y, alpha) {
  n <- nrow(X)
  q <- ncol(Y)

  # Gradient at B=0: X'Y / n
  grad <- crossprod(X, Y) / n

  # For sparse group lasso, lambda_max depends on both norms
  # Approximate: max of L2 norm of rows (group) and L_inf norm (element)

  row_l2_norms <- sqrt(rowSums(grad^2))
  max_element <- max(abs(grad))

  # Weighted combination based on alpha
  lambda_max <- max(
    max_element / alpha,           # For L1 part
    max(row_l2_norms) / (1 - alpha + 1e-10)  # For group part
  )

  return(lambda_max)
}


#' Fit Sparse Group LASSO using coordinate descent
#' @keywords internal
.fit_sgl <- function(X, Y, lambda, alpha, B_init = NULL, max_iter = 500, tol = 1e-4) {
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)

  # Initialize
  if (is.null(B_init)) {
    B <- matrix(0, nrow = p, ncol = q)
  } else {
    B <- B_init
  }

  # Precompute X'X diagonal and X'Y
  XtX_diag <- colSums(X^2) / n
  XtY <- crossprod(X, Y) / n

  # Residuals
  R <- Y - X %*% B

  lambda1 <- alpha * lambda        # L1 penalty
  lambda2 <- (1 - alpha) * lambda  # Group penalty

  for (iter in 1:max_iter) {
    B_old <- B
    max_change <- 0

    # Coordinate descent over SNPs (rows of B)
    for (j in 1:p) {
      # Current coefficient row (ensure it's a vector)
      b_j <- as.vector(B[j, ])

      # Partial residual (add back current contribution)
      # Use tcrossprod for proper n x q result: X[,j] (n x 1) times b_j (1 x q)
      R <- R + X[, j] %o% b_j

      # Gradient for this SNP (convert to vector)
      grad_j <- as.vector(crossprod(X[, j], R)) / n

      # Soft-thresholding for L1 penalty
      b_j_soft <- sign(grad_j) * pmax(abs(grad_j) - lambda1, 0)

      # Group soft-thresholding for L2 penalty
      norm_soft <- sqrt(sum(b_j_soft^2))
      if (norm_soft > 0) {
        shrink_factor <- max(0, 1 - lambda2 / (XtX_diag[j] * norm_soft + 1e-10))
        b_j_new <- as.vector(shrink_factor * b_j_soft / XtX_diag[j])
      } else {
        b_j_new <- rep(0, q)
      }

      # Update
      B[j, ] <- b_j_new
      R <- R - X[, j] %o% b_j_new

      max_change <- max(max_change, max(abs(b_j_new - b_j)))
    }

    # Check convergence
    if (max_change < tol) {
      break
    }
  }

  return(B)
}


#' Cross-validated Sparse Group LASSO with automatic alpha selection
#'
#' Runs SGL over a grid of alpha values and selects the best combination
#' of alpha and lambda via cross-validation.
#'
#' @param X matrix, design matrix of SNP dosages
#' @param Y matrix, matrix of G isoform expression across columns
#' @param alpha_seq numeric vector, sequence of alpha values to try
#' @param nlambda int, number of lambda values per alpha
#' @param nfolds int, number of CV folds
#' @param verbose logical
#' @param seed int, random seed
#' @param par logical, use parallel processing for CV folds
#' @param n.cores int, number of cores for parallel processing
#'
#' @return isotwas_model object with best alpha and lambda
#'
#' @export
multivariate_sgl_cv <- function(X,
                                 Y,
                                 alpha_seq = c(0.1, 0.5, 0.9),
                                 nlambda = 15,
                                 nfolds = 5,
                                 verbose = FALSE,
                                 seed = 123,
                                 par = FALSE,
                                 n.cores = 1) {

  best_mse <- Inf
 best_model <- NULL

  for (alpha in alpha_seq) {
    if (verbose) cat(sprintf("Trying alpha = %.2f\n", alpha))

    model <- multivariate_sgl(
      X = X, Y = Y,
      alpha = alpha,
      nlambda = nlambda,
      nfolds = nfolds,
      verbose = FALSE,
      seed = seed,
      par = par,
      n.cores = n.cores
    )

    # Compute mean MSE from CV predictions
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

  if (verbose) cat(sprintf("Best alpha: %.2f\n", best_model$alpha))

  return(best_model)
}
