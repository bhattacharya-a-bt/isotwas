#' Multivariate Super Learner Stacking
#'
#' Combines predictions from multiple base methods using optimal weights
#' learned via cross-validation. This meta-learner approach typically
#' improves upon the best single method.
#'
#' @param X matrix, design matrix of SNP dosages
#' @param Y matrix, matrix of G isoform expression across columns
#' @param Y.rep matrix, matrix of G isoform expression with replicates (optional)
#' @param R int, number of replicates (optional)
#' @param id vector, vector of sample ids (optional)
#' @param Omega matrix, precision matrix of Y (optional, computed if not provided)
#' @param base_methods character vector, methods to use as base learners.
#'   Available: "multi_enet", "univariate_enet", "ridge", "lasso", "spls"
#' @param nfolds_stack int, number of CV folds for stacking weights
#' @param verbose logical
#' @param seed int, random seed
#' @param parallel logical, use parallel processing for base methods
#' @param n.cores int, number of cores for parallel processing
#'
#' @return isotwas_model object with stacked predictions
#'
#' @export
multivariate_stacking <- function(X,
                                   Y,
                                   Y.rep = NULL,
                                   R = NULL,
                                   id = NULL,
                                   Omega = NULL,
                                   base_methods = c("multi_enet", "univariate_enet", "ridge", "lasso"),
                                   nfolds_stack = 5,
                                   verbose = FALSE,
                                   seed = 123,
                                   parallel = FALSE,
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

  # Compute Omega if needed and not provided
  if (is.null(Omega) && !is.null(Y.rep)) {
    if (verbose) cat("Computing Omega...\n")
    omega_list <- compute_omega(Y, Y.rep, R, id, method = "mean",
                                nlambda = 5, verbose = FALSE)
    Omega <- omega_list$icov[[length(omega_list$icov)]]
  } else if (is.null(Omega)) {
    Omega <- diag(q)
  }

  # Precompute optimal hyperparameters for each method using quick CV on full data
  # This avoids nested CV which is the main bottleneck
  if (verbose) cat("Pre-computing optimal hyperparameters...\n")
  param_cache <- .precompute_hyperparams(X, Y, base_methods, seed)

  # Create CV folds for stacking
  stack_folds <- caret::createFolds(1:n, k = nfolds_stack, returnTrain = TRUE)

  # Initialize prediction matrices for each base method
  base_preds <- list()
  for (method in base_methods) {
    base_preds[[method]] <- matrix(0, nrow = n, ncol = q)
  }

  if (verbose) cat("Generating base learner predictions via CV...\n")

  # Get CV predictions from each base method
  for (fold_idx in 1:nfolds_stack) {
    if (verbose) cat(sprintf("  Fold %d/%d\n", fold_idx, nfolds_stack))

    train_idx <- stack_folds[[fold_idx]]
    test_idx <- setdiff(1:n, train_idx)

    X_train <- X[train_idx, , drop = FALSE]
    Y_train <- Y[train_idx, , drop = FALSE]
    X_test <- X[test_idx, , drop = FALSE]

    # Fit all base methods (optionally in parallel)
    if (parallel && n.cores > 1 && length(base_methods) > 1) {
      fold_preds <- parallel::mclapply(base_methods, function(method) {
        .fit_base_method_fast(method, X_train, Y_train, X_test,
                              param_cache[[method]], seed + fold_idx)
      }, mc.cores = min(n.cores, length(base_methods)))
      names(fold_preds) <- base_methods
    } else {
      fold_preds <- lapply(base_methods, function(method) {
        .fit_base_method_fast(method, X_train, Y_train, X_test,
                              param_cache[[method]], seed + fold_idx)
      })
      names(fold_preds) <- base_methods
    }

    for (method in base_methods) {
      base_preds[[method]][test_idx, ] <- fold_preds[[method]]
    }
  }

  if (verbose) cat("Learning stacking weights...\n")

  # Learn stacking weights for each transcript using NNLS
  stack_weights <- matrix(0, nrow = length(base_methods), ncol = q)
  rownames(stack_weights) <- base_methods
  colnames(stack_weights) <- tx_names

  final_preds <- matrix(0, nrow = n, ncol = q)

  for (j in 1:q) {
    # Build design matrix of base predictions for this transcript
    Z <- sapply(base_methods, function(m) base_preds[[m]][, j])

    # Non-negative least squares for weights
    weights <- .nnls_weights_fast(Z, Y[, j])
    stack_weights[, j] <- weights

    # Final stacked prediction
    final_preds[, j] <- Z %*% weights
  }

  if (verbose) cat("Fitting final models on full data...\n")

  # Fit final models on full data for each base method (using cached params)
  final_base_models <- list()
  for (method in base_methods) {
    final_base_models[[method]] <- .fit_base_method_coef(
      method, X, Y, param_cache[[method]], seed
    )
  }

  # Combine base model weights using stacking weights
  stacked_weights_list <- .combine_base_weights(
    final_base_models, stack_weights, base_methods, tx_names, snp_names
  )

  # Compute R² and p-values for stacked predictions
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
    weights_df <- stacked_weights_list[[j]]

    transcripts[[tx_name]] <- create_transcript_model(
      transcript_id = tx_name,
      weights = weights_df,
      r2 = r2_vec[j],
      pvalue = p_vec[j],
      predicted = final_preds[, j]
    )
    transcripts[[tx_name]]$stack_weights <- stack_weights[, j]
  }

  result <- create_isotwas_model(
    method = "stacking",
    transcripts = transcripts,
    n_samples = n,
    n_snps = p
  )

  result$base_methods <- base_methods
  result$stack_weights_matrix <- stack_weights

  return(result)
}


#' Precompute optimal hyperparameters for each method
#' @keywords internal
.precompute_hyperparams <- function(X, Y, base_methods, seed) {
  param_cache <- list()
  q <- ncol(Y)

  for (method in base_methods) {
    set.seed(seed)

    if (method %in% c("multi_enet", "ridge", "lasso")) {
      alpha <- switch(method,
                      "multi_enet" = 0.5,
                      "ridge" = 0,
                      "lasso" = 1)

      # Quick 3-fold CV to get lambda
      fit <- glmnet::cv.glmnet(
        x = X, y = Y,
        family = "mgaussian",
        alpha = alpha,
        nfolds = 3,
        standardize = FALSE
      )
      param_cache[[method]] <- list(lambda = fit$lambda.min)

    } else if (method == "univariate_enet") {
      # For univariate, get lambda for each transcript
      lambdas <- sapply(1:q, function(j) {
        fit <- glmnet::cv.glmnet(x = X, y = Y[, j], alpha = 0.5, nfolds = 3)
        fit$lambda.min
      })
      param_cache[[method]] <- list(lambda = lambdas)

    } else if (method == "spls") {
      # Pre-compute optimal K and eta via quick CV
      # Use coarser grid for speed: K in 1:min(q,3), eta in c(0.1, 0.5, 0.9)
      max_K <- min(q, 3)
      X_clean <- X[, apply(X, 2, sd) != 0, drop = FALSE]

      opt <- tryCatch({
        spls::cv.spls(
          x = X_clean, y = Y,
          fold = 3,  # Quick 3-fold CV
          K = 1:max_K,
          eta = c(0.1, 0.5, 0.9),  # Coarser grid
          scale.x = FALSE, scale.y = FALSE,
          plot.it = FALSE
        )
      }, error = function(e) {
        list(K.opt = 1, eta.opt = 0.5)
      })

      param_cache[[method]] <- list(K = opt$K.opt, eta = opt$eta.opt)
    }
  }

  return(param_cache)
}


#' Fast base method fitting using pre-computed hyperparameters
#' @keywords internal
.fit_base_method_fast <- function(method, X_train, Y_train, X_test, params, seed) {
  q <- ncol(Y_train)
  n_test <- nrow(X_test)
  preds <- matrix(0, nrow = n_test, ncol = q)

  tryCatch({
    if (method == "multi_enet") {
      fit <- glmnet::glmnet(
        x = X_train, y = Y_train,
        family = "mgaussian",
        alpha = 0.5,
        lambda = params$lambda,
        standardize = FALSE
      )
      pred_list <- predict(fit, newx = X_test, s = params$lambda)
      for (j in 1:q) {
        preds[, j] <- pred_list[, j, 1]
      }

    } else if (method == "univariate_enet") {
      for (j in 1:q) {
        fit <- glmnet::glmnet(
          x = X_train, y = Y_train[, j],
          alpha = 0.5,
          lambda = params$lambda[j]
        )
        preds[, j] <- predict(fit, newx = X_test, s = params$lambda[j])
      }

    } else if (method == "ridge") {
      fit <- glmnet::glmnet(
        x = X_train, y = Y_train,
        family = "mgaussian",
        alpha = 0,
        lambda = params$lambda,
        standardize = FALSE
      )
      pred_list <- predict(fit, newx = X_test, s = params$lambda)
      for (j in 1:q) {
        preds[, j] <- pred_list[, j, 1]
      }

    } else if (method == "lasso") {
      fit <- glmnet::glmnet(
        x = X_train, y = Y_train,
        family = "mgaussian",
        alpha = 1,
        lambda = params$lambda,
        standardize = FALSE
      )
      pred_list <- predict(fit, newx = X_test, s = params$lambda)
      for (j in 1:q) {
        preds[, j] <- pred_list[, j, 1]
      }

    } else if (method == "spls") {
      # Filter zero-variance columns
      train_sd <- apply(X_train, 2, sd)
      test_sd <- apply(X_test, 2, sd)
      keep_cols <- which(train_sd != 0 & test_sd != 0)

      if (length(keep_cols) > 0) {
        fit <- spls::spls(
          x = X_train[, keep_cols, drop = FALSE],
          y = Y_train,
          K = params$K,
          eta = params$eta,
          scale.x = FALSE, scale.y = FALSE
        )
        preds <- predict(fit, X_test[, keep_cols, drop = FALSE])
      }
    }
  }, error = function(e) {
    # Return zeros on error
  })

  return(preds)
}


#' Get coefficients from base method using pre-computed hyperparameters
#' @keywords internal
.fit_base_method_coef <- function(method, X, Y, params, seed) {
  q <- ncol(Y)
  p <- ncol(X)
  coef_mat <- matrix(0, nrow = p, ncol = q)

  tryCatch({
    if (method == "multi_enet") {
      fit <- glmnet::glmnet(
        x = X, y = Y,
        family = "mgaussian",
        alpha = 0.5,
        lambda = params$lambda,
        standardize = FALSE
      )
      for (j in 1:q) {
        coef_mat[, j] <- as.numeric(coef(fit)[[j]])[-1]
      }

    } else if (method == "univariate_enet") {
      for (j in 1:q) {
        fit <- glmnet::glmnet(x = X, y = Y[, j], alpha = 0.5, lambda = params$lambda[j])
        coef_mat[, j] <- as.numeric(coef(fit))[-1]
      }

    } else if (method == "ridge") {
      fit <- glmnet::glmnet(
        x = X, y = Y,
        family = "mgaussian",
        alpha = 0,
        lambda = params$lambda,
        standardize = FALSE
      )
      for (j in 1:q) {
        coef_mat[, j] <- as.numeric(coef(fit)[[j]])[-1]
      }

    } else if (method == "lasso") {
      fit <- glmnet::glmnet(
        x = X, y = Y,
        family = "mgaussian",
        alpha = 1,
        lambda = params$lambda,
        standardize = FALSE
      )
      for (j in 1:q) {
        coef_mat[, j] <- as.numeric(coef(fit)[[j]])[-1]
      }

    } else if (method == "spls") {
      # Filter zero-variance columns
      keep_cols <- which(apply(X, 2, sd) != 0)

      if (length(keep_cols) > 0) {
        fit <- spls::spls(
          x = X[, keep_cols, drop = FALSE],
          y = Y,
          K = params$K,
          eta = params$eta,
          scale.x = FALSE, scale.y = FALSE
        )
        coef_spls <- coef(fit)
        coef_mat[keep_cols, ] <- coef_spls
      }
    }
  }, error = function(e) {
    # Return zeros on error
  })

  return(coef_mat)
}


#' Combine base model weights using stacking weights
#' @keywords internal
.combine_base_weights <- function(final_base_models, stack_weights,
                                   base_methods, tx_names, snp_names) {
  q <- length(tx_names)
  p <- length(snp_names)

  result <- list()

  for (j in 1:q) {
    combined_coef <- rep(0, p)

    for (m_idx in seq_along(base_methods)) {
      method <- base_methods[m_idx]
      w <- stack_weights[m_idx, j]
      if (w > 0) {
        combined_coef <- combined_coef + w * final_base_models[[method]][, j]
      }
    }

    weights_df <- tibble::tibble(
      SNP = snp_names,
      Weight = combined_coef
    )
    weights_df <- subset(weights_df, Weight != 0)

    result[[j]] <- weights_df
  }

  return(result)
}


#' Fast non-negative least squares for stacking weights
#' @keywords internal
.nnls_weights_fast <- function(Z, y) {
  n_methods <- ncol(Z)

  # Simple closed-form solution with projection to non-negative
  # Faster than iterative optimization for small n_methods
  ZtZ <- crossprod(Z)
  Zty <- crossprod(Z, y)

  # Try direct solve first
  weights <- tryCatch({
    solve(ZtZ + diag(1e-6, n_methods), Zty)
  }, error = function(e) {
    rep(1 / n_methods, n_methods)
  })

  # Project to non-negative
  weights <- pmax(weights, 0)

  # Normalize to sum to 1
  if (sum(weights) > 0) {
    weights <- weights / sum(weights)
  } else {
    weights <- rep(1 / n_methods, n_methods)
  }

  return(as.vector(weights))
}
