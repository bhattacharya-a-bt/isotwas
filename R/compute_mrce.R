#' Multivariate regression with covariance estimate
#'
#' The function trains multivariate lasso with a given covariance estimate
#'
#' @param X matrix, design matrix of SNP dosages
#' @param Y matrix, matrix of G isoform expression across columns
#' @param lambda vector, lambda penalty vector for LASSO to tune on
#' @param Omega matrix, precision matrix of Y
#' @param nfolds int, number of CV folds
#' @param tol.in numeric, tolerance for objective difference
#' @param maxit.in int, maximum number of iteractions
#' @param verbose logical
#' @param seed int, random seed
#' @param nlambda int, number of lambda value for parameter selection
#' @param parallel logical, whether to parallelize CV folds
#' @param n.cores int, number of cores for parallelization
#'
#' @return CV MRCE fit
#'
#' @export
compute_mrce = function(X,
                        Y,
                        lambda = NULL,
                        nlambda = 20,
                        Omega,
                        nfolds = 5,
                        tol.in,
                        maxit.in,
                        verbose,
                        seed,
                        parallel = FALSE,
                        n.cores = NULL){

  n = nrow(X)
  p = ncol(X)
  q = ncol(Y)

  # Smarter lambda grid: data-driven range

  if (is.null(lambda)){
    # Compute lambda_max (smallest lambda that gives all zeros)
    mx = colMeans(X)
    my = colMeans(Y)
    X.centered = scale(X, center = mx, scale = FALSE)
    Y.centered = scale(Y, center = my, scale = FALSE)
    xty = crossprod(X.centered, Y.centered)
    xtyom = xty %*% Omega
    lambda_max = max(abs(xtyom)) / n
    lambda_min = lambda_max * 1e-4
    lambda = 10^seq(log10(lambda_max), log10(lambda_min), length.out = nlambda)
  }

  set.seed(seed)
  train.folds = caret::createFolds(1:nrow(Y),
                                   k = nfolds,
                                   returnTrain = TRUE)

  # Precompute fold-specific data (variance filtering, matrices)
  fold_data = lapply(1:nfolds, function(tr) {
    train_idx = train.folds[[tr]]
    test_idx = setdiff(1:n, train_idx)

    Y.tr = Y[train_idx, , drop = FALSE]
    X.tr = X[train_idx, , drop = FALSE]
    X.test = X[test_idx, , drop = FALSE]

    # Variance filter - compute once per fold
    sd_train = apply(X.tr, 2, sd)
    sd_test = apply(X.test, 2, sd)
    valid_snps = which(sd_train != 0 & sd_test != 0)

    X.tr = X.tr[, valid_snps, drop = FALSE]
    X.test = X.test[, valid_snps, drop = FALSE]

    # Precompute matrices for this fold
    mx = colMeans(X.tr)
    my = colMeans(Y.tr)
    X.tr.centered = scale(X.tr, center = mx, scale = FALSE)
    Y.tr.centered = scale(Y.tr, center = my, scale = FALSE)

    xtx = crossprod(X.tr.centered)
    xty = crossprod(X.tr.centered, Y.tr.centered)
    xtyom = xty %*% Omega
    tolmult = sum(crossprod(Y.tr.centered) * Omega)

    list(
      train_idx = train_idx,
      test_idx = test_idx,
      valid_snps = valid_snps,
      X.test = X.test,
      n_train = nrow(X.tr),
      mx = mx,
      my = my,
      xtx = xtx,
      xtyom = xtyom,
      tolmult = tolmult
    )
  })

  # Function to run CV for a single fold with warm starts across lambda path
  run_fold_cv = function(fold_idx, fold_data, lambda, Omega, tol.in, maxit.in) {
    fd = fold_data[[fold_idx]]
    n_lambda = length(lambda)
    p_fold = length(fd$valid_snps)

    # Initialize coefficient matrix for warm starts
    B_warm = matrix(0, nrow = p_fold, ncol = q)

    # Store predictions for each lambda
    fold_preds = array(0, dim = c(length(fd$test_idx), q, n_lambda))

    for (lam_idx in 1:n_lambda) {
      # Use warm start from previous lambda
      mod = compute_fixed_precomputed(
        xtx = fd$xtx,
        xtyom = fd$xtyom,
        tolmult = fd$tolmult,
        n = fd$n_train,
        p = p_fold,
        q = q,
        lam2 = lambda[lam_idx],
        Omega = Omega,
        tol.in = tol.in,
        maxit.in = maxit.in,
        silent = !verbose,
        B0 = B_warm,
        mx = fd$mx,
        my = fd$my
      )

      B_warm = mod$Bhat  # Warm start for next lambda
      fold_preds[, , lam_idx] = fd$X.test %*% mod$Bhat
    }

    list(test_idx = fd$test_idx, preds = fold_preds)
  }

  # Run CV - either parallel or sequential
  if (parallel && !is.null(n.cores) && n.cores > 1) {
    fold_results = parallel::mclapply(
      1:nfolds,
      function(tr) run_fold_cv(tr, fold_data, lambda, Omega, tol.in, maxit.in),
      mc.cores = min(n.cores, nfolds)
    )
  } else {
    fold_results = lapply(
      1:nfolds,
      function(tr) {
        if (verbose) cat("Fold", tr, "of", nfolds, "\n")
        run_fold_cv(tr, fold_data, lambda, Omega, tol.in, maxit.in)
      }
    )
  }

  # Assemble predictions from all folds
  pred_array = array(0, dim = c(n, q, length(lambda)))
  for (fr in fold_results) {
    pred_array[fr$test_idx, , ] = fr$preds
  }

  # Vectorized R² and p-value calculation (much faster than lm())
  r2Mat = matrix(nrow = q, ncol = length(lambda))
  pMat = matrix(nrow = q, ncol = length(lambda))

  for (lam_idx in 1:length(lambda)) {
    pred = pred_array[, , lam_idx]
    pred[is.na(pred)] = 0

    for (tx in 1:q) {
      y_obs = Y[, tx]
      y_pred = pred[, tx]

      # Fast R² calculation
      ss_res = sum((y_obs - y_pred)^2)
      ss_tot = sum((y_obs - mean(y_obs))^2)

      if (ss_tot > 0) {
        r2_raw = 1 - ss_res / ss_tot
        # Adjusted R²
        r2_adj = 1 - (1 - r2_raw) * (n - 1) / (n - 2)
        r2Mat[tx, lam_idx] = max(0, r2_adj)
      } else {
        r2Mat[tx, lam_idx] = 0
      }

      # Fast p-value from correlation
      if (sd(y_pred) > 0) {
        r = cor(y_obs, y_pred)
        t_stat = r * sqrt((n - 2) / (1 - r^2))
        pMat[tx, lam_idx] = 2 * pt(-abs(t_stat), df = n - 2)
      } else {
        pMat[tx, lam_idx] = 1
      }
    }
  }

  # Select best lambda
  mean_r2 = colMeans(r2Mat, na.rm = TRUE)
  best_lambda_idx = which.max(mean_r2)
  best_lambda = lambda[best_lambda_idx]

  # Get final predictions at best lambda
  pred = pred_array[, , best_lambda_idx]
  pred[is.na(pred)] = 0

  # Fit final model on all data with warm start
  mx_full = colMeans(X)
  my_full = colMeans(Y)
  X.centered = scale(X, center = mx_full, scale = FALSE)
  Y.centered = scale(Y, center = my_full, scale = FALSE)
  xtx_full = crossprod(X.centered)
  xty_full = crossprod(X.centered, Y.centered)
  xtyom_full = xty_full %*% Omega
  tolmult_full = sum(crossprod(Y.centered) * Omega)

  final.model = compute_fixed_precomputed(
    xtx = xtx_full,
    xtyom = xtyom_full,
    tolmult = tolmult_full,
    n = n,
    p = p,
    q = q,
    lam2 = best_lambda,
    Omega = Omega,
    tol.in = tol.in,
    maxit.in = maxit.in,
    silent = !verbose,
    B0 = NULL,
    mx = mx_full,
    my = my_full
  )

  # Build output using isotwas_model class
  transcripts <- list()
  for (i in 1:q) {
    tx_name <- colnames(Y)[i]
    weights <- tibble::tibble(
      SNP = colnames(X),
      Weight = final.model$Bhat[, i]
    )
    weights <- subset(weights, Weight != 0)

    transcripts[[tx_name]] <- create_transcript_model(
      transcript_id = tx_name,
      weights = weights,
      r2 = r2Mat[i, best_lambda_idx],
      pvalue = pMat[i, best_lambda_idx],
      predicted = pred[, i]
    )
  }

  result <- create_isotwas_model(
    method = "mrce_lasso",
    transcripts = transcripts,
    n_samples = n,
    n_snps = p
  )

  return(result)
}
