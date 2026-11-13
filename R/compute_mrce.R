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
#' @param r2_method character, how the cross-validated R-squared is scored:
#'   "refit" (the default), the adjusted R-squared of lm(obs ~ pred) that
#'   calc.r2() and the other methods use, or "oos", the direct out-of-sample
#'   R-squared 1 - SS_res/SS_tot that this method used before.
#'   compute_isotwas() passes "oos" under \code{r2_scale = 'legacy'}.
#'   The mean of this score across transcripts also chooses lambda.
#' @param floor_r2 logical, clamp a negative R-squared to 0, as earlier
#'   versions did. Default FALSE.
#' @param fold_cache list or NULL, the Omega-independent quantities from an
#'   internal helper, so that fits to several Omegas on the same X, Y, nfolds
#'   and seed compute them once. Leave NULL.
#' @param folds optional cross-validation folds: a list with one vector of
#'   training-row indices per fold, or a vector of fold labels 1..K, one per
#'   sample. NULL (the default) draws them from \code{seed} as before;
#'   compute_isotwas passes every method the same folds.
#'
#' @return an isotwas_model, with the lambda path in \code{lambda}, the chosen
#'   value in \code{best_lambda}, and in \code{cv_r2} the cross-validated
#'   R-squared of each transcript (rows) at each lambda (columns), whose column
#'   means chose it
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
                        n.cores = NULL,
                        r2_method = c("refit", "oos"),
                        floor_r2 = FALSE,
                        fold_cache = NULL,
                        folds = NULL){

  r2_method = match.arg(r2_method)

  n = nrow(X)
  p = ncol(X)
  q = ncol(Y)

  # Smarter lambda grid: data-driven range

  # Everything that does not depend on Omega: the folds, each fold's
  # variance filter and cross-products, and the full-data cross-products.
  # compute_isotwas fits one Omega after another on the same data, so it
  # builds this once and passes it in.
  if (is.null(fold_cache)) {
    fold_cache = .mrce_fold_cache(X, Y, nfolds, seed, folds)
  }
  set.seed(seed)
  train.folds = caret::createFolds(1:nrow(Y),
                                   k = nfolds,
                                   returnTrain = TRUE)
  if (!is.null(folds)) {
    folds = .as_folds(folds, nrow(Y))
    train.folds = folds$train
    nfolds = folds$nfolds
  }
  if (!identical(train.folds, fold_cache$train.folds)) {
    stop("fold_cache was built with different folds")
  }
  full = fold_cache$full

  if (is.null(lambda)){
    # Compute lambda_max (smallest lambda that gives all zeros)
    xtyom = full$xty %*% Omega
    lambda_max = max(abs(xtyom)) / n
    lambda_min = lambda_max * 1e-4
    lambda = 10^seq(log10(lambda_max), log10(lambda_min), length.out = nlambda)
  }

  fold_data = lapply(fold_cache$folds, function(fd) {
    fd$xtyom = fd$xty %*% Omega
    fd$tolmult = sum(fd$yty * Omega)
    fd
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
      # SNPs with no weight add exact zeros to the prediction, so skip them.
      # The fit is on centred data, so its intercept muhat = my - mx'B goes
      # back in
      nz = which(rowSums(mod$Bhat != 0) > 0)
      fold_preds[, , lam_idx] = fd$X.test[, nz, drop = FALSE] %*%
        mod$Bhat[nz, , drop = FALSE] +
        rep(mod$muhat, each = length(fd$test_idx))
    }

    list(test_idx = fd$test_idx, preds = fold_preds)
  }

  # Run CV - either parallel or sequential
  if (parallel && !is.null(n.cores) && n.cores > 1) {
    fold_results = .mclapply_seeded(
      1:nfolds,
      function(tr) run_fold_cv(tr, fold_data, lambda, Omega, tol.in, maxit.in),
      mc.cores = min(n.cores, nfolds),
      seed = seed
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

      # P is needed only at the selected lambda, scored below
      scored = .score_r2(y_obs, y_pred, method = r2_method, pvalue = FALSE)
      r2Mat[tx, lam_idx] = if (floor_r2) max(0, scored$R2) else scored$R2
    }
  }

  # Select best lambda
  mean_r2 = colMeans(r2Mat, na.rm = TRUE)
  best_lambda_idx = which.max(mean_r2)
  best_lambda = lambda[best_lambda_idx]

  # Get final predictions at best lambda
  pred = pred_array[, , best_lambda_idx]
  pred[is.na(pred)] = 0
  for (tx in 1:q) {
    pMat[tx, best_lambda_idx] = .score_r2(Y[, tx], pred[, tx],
                                          method = r2_method)$P
  }

  # Fit final model on all data with warm start
  final.model = compute_fixed_precomputed(
    xtx = full$xtx,
    xtyom = full$xty %*% Omega,
    tolmult = sum(full$yty * Omega),
    n = n,
    p = p,
    q = q,
    lam2 = best_lambda,
    Omega = Omega,
    tol.in = tol.in,
    maxit.in = maxit.in,
    silent = !verbose,
    B0 = NULL,
    mx = full$mx,
    my = full$my
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
  # The lambda path and each transcript's cross-validated R2 along it, from
  # which best_lambda was chosen
  result$lambda <- lambda
  result$best_lambda <- best_lambda
  result$cv_r2 <- r2Mat

  result <- .drop_null_fits(result, Y)

  return(result)
}

#' Omega-independent part of compute_mrce
#'
#' The folds compute_mrce draws for this seed, each fold's variance filter,
#' centring and cross-products, and the full-data cross-products. The RNG
#' state is restored on exit, since compute_mrce draws the folds again itself.
#'
#' @return list(train.folds, folds, full)
#' @noRd
.mrce_fold_cache <- function(X, Y, nfolds, seed, folds = NULL) {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  }
  n <- nrow(X)
  set.seed(seed)
  train.folds <- caret::createFolds(1:nrow(Y), k = nfolds, returnTrain = TRUE)
  if (!is.null(folds)) train.folds <- .as_folds(folds, nrow(Y))$train

  folds <- lapply(seq_along(train.folds), function(tr) {
    train_idx <- train.folds[[tr]]
    test_idx <- setdiff(1:n, train_idx)

    Y.tr <- Y[train_idx, , drop = FALSE]
    X.tr <- X[train_idx, , drop = FALSE]
    X.test <- X[test_idx, , drop = FALSE]

    # Variance filter
    sd_train <- apply(X.tr, 2, sd)
    sd_test <- apply(X.test, 2, sd)
    valid_snps <- which(sd_train != 0 & sd_test != 0)

    X.tr <- X.tr[, valid_snps, drop = FALSE]
    X.test <- X.test[, valid_snps, drop = FALSE]

    mx <- colMeans(X.tr)
    my <- colMeans(Y.tr)
    X.tr.centered <- scale(X.tr, center = mx, scale = FALSE)
    Y.tr.centered <- scale(Y.tr, center = my, scale = FALSE)

    list(train_idx = train_idx,
         test_idx = test_idx,
         valid_snps = valid_snps,
         X.test = X.test,
         n_train = nrow(X.tr),
         mx = mx,
         my = my,
         xtx = crossprod(X.tr.centered),
         xty = crossprod(X.tr.centered, Y.tr.centered),
         yty = crossprod(Y.tr.centered))
  })

  mx <- colMeans(X)
  my <- colMeans(Y)
  X.centered <- scale(X, center = mx, scale = FALSE)
  Y.centered <- scale(Y, center = my, scale = FALSE)

  list(train.folds = train.folds,
       folds = folds,
       full = list(mx = mx,
                   my = my,
                   xtx = crossprod(X.centered),
                   xty = crossprod(X.centered, Y.centered),
                   yty = crossprod(Y.centered)))
}
