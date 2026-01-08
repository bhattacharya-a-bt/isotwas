#' Multivariate Graph-Regularized Regression
#'
#' Fits a multivariate regression model with graph Laplacian regularization.
#' This method encourages isoforms that are similar (e.g., share exons) to have
#' similar prediction weight vectors.
#'
#' The optimization problem is:
#' \deqn{\min_B \frac{1}{2n}||Y - XB||_F^2 + \lambda_1 ||B||_1 + \frac{\lambda_g}{2} \text{tr}(B^T B L)}
#'
#' where L is the graph Laplacian computed from the similarity matrix:
#' \deqn{L = D - S}
#' with D being the degree matrix (diagonal with row sums of S).
#'
#' The graph penalty encourages isoforms with high similarity to have similar
#' coefficient vectors across SNPs.
#'
#' @param X matrix, design matrix of SNP dosages (n x p)
#' @param Y matrix, matrix of isoform expression across columns (n x q)
#' @param similarity_matrix matrix, q x q symmetric matrix of pairwise isoform
#'   similarities (e.g., Jaccard index of shared exons). Values should be >= 0.
#'   If NULL, uses identity (no graph regularization, equivalent to lasso).
#' @param lambda1 numeric or NULL, L1 penalty parameter. If NULL, selected by CV.
#' @param lambda_graph numeric or NULL, graph penalty parameter. If NULL, selected by CV.
#' @param alpha numeric, elastic net mixing (0 = ridge, 1 = lasso). Default 0.5.
#' @param nlambda int, number of lambda1 values to try in CV
#' @param nlambda_graph int, number of lambda_graph values to try in CV
#' @param lambda_graph_seq numeric vector, specific lambda_graph values to try.
#'   If NULL, auto-generated.
#' @param normalize_laplacian logical, use normalized Laplacian. Default FALSE.
#' @param nfolds int, number of CV folds
#' @param standardize logical, standardize X before fitting
#' @param verbose logical, print progress
#' @param seed int, random seed
#'
#' @return isotwas_model object containing:
#'   \itemize{
#'     \item transcripts: list of transcript_model objects
#'     \item best_lambda1: optimal L1 penalty
#'     \item best_lambda_graph: optimal graph penalty
#'     \item laplacian: the Laplacian matrix used
#'   }
#'
#' @details
#' The similarity matrix encodes prior knowledge about isoform relationships.
#' For isoform-level TWAS, this is typically derived from shared exon structure:
#' isoforms sharing more exons are expected to have more similar cis-regulatory
#' effects, as they share more of the same genetic signal.
#'
#' @export
multivariate_graph_reg <- function(X,
                                    Y,
                                    similarity_matrix = NULL,
                                    lambda1 = NULL,
                                    lambda_graph = NULL,
                                    alpha = 0.5,
                                    nlambda = 15,
                                    nlambda_graph = 5,
                                    lambda_graph_seq = NULL,
                                    normalize_laplacian = FALSE,
                                    nfolds = 5,
                                    standardize = FALSE,
                                    verbose = FALSE,
                                    seed = 123) {

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

  # Compute Laplacian from similarity matrix
  if (is.null(similarity_matrix)) {
    # No graph regularization - use identity (makes graph penalty = ridge on B)
    L <- diag(q)
    if (verbose) cat("No similarity matrix provided - using identity (no graph structure)\n")
  } else {
    # Validate similarity matrix
    if (!is.matrix(similarity_matrix)) {
      similarity_matrix <- as.matrix(similarity_matrix)
    }
    if (nrow(similarity_matrix) != q || ncol(similarity_matrix) != q) {
      stop("similarity_matrix must be q x q where q = number of isoforms (", q, ")")
    }
    if (!isSymmetric(similarity_matrix, tol = 1e-8)) {
      warning("similarity_matrix is not symmetric - symmetrizing")
      similarity_matrix <- (similarity_matrix + t(similarity_matrix)) / 2
    }
    if (any(similarity_matrix < 0)) {
      warning("similarity_matrix contains negative values - taking absolute value")
      similarity_matrix <- abs(similarity_matrix)
    }

    # Compute Laplacian: L = D - S
    L <- .compute_laplacian(similarity_matrix, normalize = normalize_laplacian)
    if (verbose) cat("Computed Laplacian from similarity matrix\n")
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

  # Precompute X'X and X'Y
  XtX <- crossprod(X_scaled) / n
  XtY <- crossprod(X_scaled, Y_centered) / n

  # Compute lambda1 sequence (L1 penalty)
  lambda1_max <- max(abs(XtY)) / alpha
  lambda1_min <- lambda1_max * 0.01
  lambda1_seq <- exp(seq(log(lambda1_max), log(lambda1_min), length.out = nlambda))

  # Compute lambda_graph sequence
  if (is.null(lambda_graph_seq)) {
    # Scale based on Laplacian eigenvalues
    L_eigs <- eigen(L, symmetric = TRUE, only.values = TRUE)$values
    max_eig <- max(L_eigs)
    if (max_eig > 0) {
      lambda_graph_seq <- c(0, 10^seq(-2, 1, length.out = nlambda_graph - 1) * max_eig)
    } else {
      lambda_graph_seq <- c(0)
    }
  }

  # Create CV folds
  cv_folds <- caret::createFolds(1:n, k = nfolds, returnTrain = TRUE)

  # If lambdas provided, skip CV
  if (!is.null(lambda1) && !is.null(lambda_graph)) {
    best_lambda1 <- lambda1
    best_lambda_graph <- lambda_graph
    if (verbose) cat(sprintf("Using provided lambda1=%.4f, lambda_graph=%.4f\n",
                             best_lambda1, best_lambda_graph))
  } else {
    # Cross-validation over lambda1 and lambda_graph grid
    if (verbose) cat("Running cross-validation for parameter selection...\n")

    cv_results <- matrix(Inf, nrow = length(lambda1_seq), ncol = length(lambda_graph_seq))

    for (lg_idx in seq_along(lambda_graph_seq)) {
      lg <- lambda_graph_seq[lg_idx]
      if (verbose) cat(sprintf("  lambda_graph = %.4f (%d/%d)\n",
                               lg, lg_idx, length(lambda_graph_seq)))

      for (fold_idx in 1:nfolds) {
        train_idx <- cv_folds[[fold_idx]]
        test_idx <- setdiff(1:n, train_idx)

        X_train <- X_scaled[train_idx, , drop = FALSE]
        Y_train <- Y_centered[train_idx, , drop = FALSE]
        X_test <- X_scaled[test_idx, , drop = FALSE]
        Y_test <- Y_centered[test_idx, , drop = FALSE]

        # Warm start path
        B_warm <- matrix(0, p, q)

        for (l1_idx in seq_along(lambda1_seq)) {
          l1 <- lambda1_seq[l1_idx]

          B_warm <- .fit_graph_reg(X_train, Y_train, L, l1, lg, alpha,
                                   B_init = B_warm, max_iter = 200, tol = 1e-4)

          pred <- X_test %*% B_warm
          mse <- mean((Y_test - pred)^2)
          cv_results[l1_idx, lg_idx] <- cv_results[l1_idx, lg_idx] + mse / nfolds
        }
      }
    }

    # Find best parameters
    best_idx <- which(cv_results == min(cv_results), arr.ind = TRUE)[1, ]
    best_lambda1 <- lambda1_seq[best_idx[1]]
    best_lambda_graph <- lambda_graph_seq[best_idx[2]]

    if (verbose) {
      cat(sprintf("Best lambda1: %.6f, Best lambda_graph: %.6f\n",
                  best_lambda1, best_lambda_graph))
    }
  }

  # Fit final model
  if (verbose) cat("Fitting final model...\n")
  B_final <- .fit_graph_reg(X_scaled, Y_centered, L, best_lambda1, best_lambda_graph,
                            alpha, B_init = NULL, max_iter = 1000, tol = 1e-5)

  # Rescale coefficients if standardized
  if (standardize) {
    B_final <- B_final / X_sds
  }

  # Get CV predictions for R² calculation
  cv_preds <- matrix(0, nrow = n, ncol = q)
  for (fold_idx in 1:nfolds) {
    train_idx <- cv_folds[[fold_idx]]
    test_idx <- setdiff(1:n, train_idx)

    X_train <- X_scaled[train_idx, , drop = FALSE]
    Y_train <- Y_centered[train_idx, , drop = FALSE]
    X_test <- X_scaled[test_idx, , drop = FALSE]

    B_fold <- .fit_graph_reg(X_train, Y_train, L, best_lambda1, best_lambda_graph,
                             alpha, B_init = NULL, max_iter = 500, tol = 1e-4)
    cv_preds[test_idx, ] <- X_test %*% B_fold
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
    method = "graph_regularized",
    transcripts = transcripts,
    n_samples = n,
    n_snps = p
  )

  # Add graph-reg specific info
  result$best_lambda1 <- best_lambda1
  result$best_lambda_graph <- best_lambda_graph
  result$alpha <- alpha
  result$laplacian <- L
  result$similarity_matrix <- similarity_matrix
  if (exists("cv_results")) {
    result$cv_results <- cv_results
    result$lambda1_seq <- lambda1_seq
    result$lambda_graph_seq <- lambda_graph_seq
  }

  return(result)
}


#' Compute graph Laplacian from similarity matrix
#'
#' @param S similarity matrix (symmetric, non-negative)
#' @param normalize logical, compute normalized Laplacian
#' @return Laplacian matrix
#' @keywords internal
.compute_laplacian <- function(S, normalize = FALSE) {
  # Degree matrix
  d <- rowSums(S)

  if (normalize) {
    # Normalized Laplacian: L = I - D^{-1/2} S D^{-1/2}
    d_inv_sqrt <- ifelse(d > 0, 1 / sqrt(d), 0)
    D_inv_sqrt <- diag(d_inv_sqrt)
    L <- diag(nrow(S)) - D_inv_sqrt %*% S %*% D_inv_sqrt
  } else {
    # Unnormalized Laplacian: L = D - S
    L <- diag(d) - S
  }

  return(L)
}


#' Fit graph-regularized elastic net using coordinate descent
#'
#' Solves: min (1/2n)||Y - XB||²_F + λ₁α||B||₁ + λ₁(1-α)/2||B||²_F + λ_g/2 tr(B'BL)
#'
#' @keywords internal
.fit_graph_reg <- function(X, Y, L, lambda1, lambda_graph, alpha,
                           B_init = NULL, max_iter = 500, tol = 1e-4) {
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)

  # Initialize
  if (is.null(B_init)) {
    B <- matrix(0, p, q)
  } else {
    B <- B_init
  }

  # Precompute
  XtX_diag <- colSums(X^2) / n
  XtY <- crossprod(X, Y) / n

  # Penalty parameters
  l1_pen <- lambda1 * alpha          # L1 penalty
  l2_pen <- lambda1 * (1 - alpha)    # Ridge penalty

  # Residual
  R <- Y - X %*% B

  for (iter in 1:max_iter) {
    B_old <- B
    max_change <- 0

    # Coordinate descent over SNPs (rows of B)
    for (j in 1:p) {
      b_j <- B[j, ]

      # Add back contribution of current SNP
      R <- R + outer(X[, j], b_j)

      # Gradient: (1/n) X_j' R
      grad_j <- as.vector(crossprod(X[, j], R)) / n

      # For graph regularization, the update involves L
      # The penalty gradient for tr(B'BL) w.r.t. b_j is 2 * B_j L
      # where B_j is the j-th row of B

      # Solve: b_j = argmin (1/2)*XtX_jj*||b||² - grad_j'b + l1*||b||₁ + (l2/2)||b||² + (lg/2)*b'Lb
      # This is: (XtX_jj + l2)*I + lg*L  as the quadratic term
      # With gradient term and L1

      if (lambda_graph > 0 && max(abs(L - diag(diag(L)))) > 1e-10) {
        # Non-trivial graph structure - need to solve coupled system
        # (XtX_jj + l2)*b + lg*L*b = grad_j, then soft-threshold

        # Form the system matrix
        H <- (XtX_diag[j] + l2_pen) * diag(q) + lambda_graph * L

        # Solve for the unconstrained optimum
        b_unconstrained <- solve(H, grad_j)

        # Apply soft-thresholding for L1
        b_j_new <- sign(b_unconstrained) * pmax(abs(b_unconstrained) - l1_pen / (XtX_diag[j] + l2_pen + lambda_graph), 0)
      } else {
        # No graph structure or identity - standard elastic net update
        denom <- XtX_diag[j] + l2_pen

        # Soft-thresholding
        b_j_new <- sign(grad_j) * pmax(abs(grad_j) - l1_pen, 0) / denom
      }

      # Update
      B[j, ] <- b_j_new
      R <- R - outer(X[, j], b_j_new)

      max_change <- max(max_change, max(abs(b_j_new - b_j)))
    }

    if (max_change < tol) {
      break
    }
  }

  return(B)
}


#' Create similarity matrix from shared exons
#'
#' Helper function to create an isoform similarity matrix based on shared exons.
#' Uses Jaccard similarity: |A ∩ B| / |A ∪ B|
#'
#' @param exon_list named list where each element is a character vector of exon IDs
#'   for that isoform. Names should match isoform/transcript IDs.
#' @param method character, similarity method: "jaccard" (default), "overlap", or "count"
#'
#' @return symmetric similarity matrix
#'
#' @details
#' Methods:
#' \itemize{
#'   \item "jaccard": Jaccard index = |shared| / |union|
#'   \item "overlap": Overlap coefficient = |shared| / min(|A|, |B|)
#'   \item "count": Raw count of shared exons
#' }
#'
#' @export
create_similarity_from_exons <- function(exon_list, method = c("jaccard", "overlap", "count")) {
  method <- match.arg(method)

  n <- length(exon_list)
  tx_names <- names(exon_list)
  if (is.null(tx_names)) {
    tx_names <- paste0("Transcript", 1:n)
  }

  S <- matrix(0, n, n)
  rownames(S) <- colnames(S) <- tx_names

  for (i in 1:n) {
    for (j in i:n) {
      exons_i <- exon_list[[i]]
      exons_j <- exon_list[[j]]

      shared <- length(intersect(exons_i, exons_j))

      sim <- switch(method,
        jaccard = shared / length(union(exons_i, exons_j)),
        overlap = shared / min(length(exons_i), length(exons_j)),
        count = shared
      )

      S[i, j] <- sim
      S[j, i] <- sim
    }
  }

  return(S)
}
