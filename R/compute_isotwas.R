#' Compute isoTWAS model for a set of isoforms/transcripts
#'
#' The function runs multivariate regression models for a set of transcripts
#' and outputs the best model based on cross-validated R-squared.
#'
#' @param X matrix, design matrix of SNP dosages
#' @param Y matrix, matrix of G isoform expression across columns
#' @param Y.rep matrix, matrix of G isoform expression with replicates
#' @param R int, number of replicates
#' @param id vector, vector of sample ids showing rep to id
#' @param omega_est character, passed to compute_omega() as its \code{method}:
#'   'replicates', 'residuals' or 'mean'
#' @param r2_scale character, how methods are scored for comparison:
#'   'coherent' (the default) or 'legacy'. See Details.
#' @param omega_nlambda int, number of omegas to generate
#' @param method character, vector of methods to use. Available methods:
#'   \itemize{
#'     \item "mrce_lasso": MRCE with lasso penalty
#'     \item "multi_enet": Multivariate elastic net
#'     \item "joinet": Joinet stacked elastic net
#'     \item "spls": Sparse partial least squares
#'     \item "sgl": Sparse group lasso
#'     \item "mtlasso": Multi-task lasso (L21 regularization)
#'     \item "stacking": Super learner stacking of multiple methods
#'     \item "graph_reg": Graph-regularized regression (requires similarity_matrix)
#'     \item "univariate": Best of univariate elastic net, BLUP, SuSiE
#'   }
#' @param similarity_matrix matrix, optional similarity matrix for graph_reg method.
#'   Can be created using \code{\link{similarity_from_gtf}} or
#'   \code{\link{create_similarity_from_exons}}. If NULL, graph_reg is skipped.
#' @param predict_nlambda int, number of lambdas in MRCE
#' @param family character, glmnet family
#' @param scale logical, T/F to scale Y by Omega
#' @param alpha numeric, elastic net mixing parameter
#' @param nfolds int, number of CV folds
#' @param verbose logical, print progress messages
#' @param par logical, uses mclapply to parallelize model fit
#' @param n.cores int, number of parallel cores
#' @param tx_names vector, character vector of tx names - order of columns of Y
#' @param seed int, random seed
#' @param return_all logical, return R2 for all models?
#' @param tol.in numeric, tolerance for objective difference
#' @param maxit.in int, maximum number of iterations
#' @param joinet_nfolds_int int, folds for joinet's inner cross-validation.
#' @param sgl_nlambda int, lambda path length for multivariate_sgl
#' @param gene_exp vector, vector of total gene expression
#' @param run_all logical, run all methods
#' @param shared_folds logical, give every method the same cross-validation
#'   folds (default TRUE), so that the cross-validated R2 compared across
#'   methods come from the same splits. The folds are those drawn by
#'   set.seed(seed) and caret::createFolds(1:N, k = nfolds,
#'   returnTrain = TRUE). FALSE lets each method draw its own, as before.
#' @param folds optional folds to use instead: a list with one vector of
#'   training-row indices per fold, or a vector of fold labels 1..K, one per
#'   sample. They override nfolds and shared_folds.
#'
#' @importFrom glmnet glmnet
#' @importFrom glmnet cv.glmnet
#' @importFrom rlist list.append
#' @importFrom caret createFolds
#' @importFrom stats lm
#' @importFrom stats coef
#'
#' @details
#' The best model per transcript is the which.max of cross-validated R-squared
#' across methods. Under the default \code{r2_scale = 'coherent'}, every
#' method scores its cross-validated prediction the same way, by the adjusted
#' R-squared of \code{lm(obs ~ pred)} (\code{r2_method = 'refit'} in the
#' fitting functions). A prediction that is constant or not positively
#' correlated with expression scores 0 and its model is emptied. A model that
#' selected no SNPs scores -Inf, and a transcript for which every model is
#' empty is labelled NA.
#'
#' \code{r2_scale = 'legacy'} reproduces the scoring of earlier versions.
#' mrce_lasso, sgl, mtlasso, stacking and graph_reg score the out-of-sample
#' R-squared 1 - SS_res/SS_tot (\code{r2_method = 'oos'}) and clamp it at 0,
#' while the other methods are scored by \code{lm(obs ~ pred)} and can score
#' below 0, and an empty model scores 0. On the same prediction the
#' out-of-sample R-squared is never higher, so selection favours the other
#' methods, and on a transcript no method predicts a clamped or empty model
#' wins.
#'
#' joinet is the most expensive method at moderate transcript counts, and most
#' of that time is its inner cross-validation. Lowering
#' \code{joinet_nfolds_int} is the first way to make the method set cheaper.
#'
#' With shared folds (the default) or folds given in \code{folds}, the
#' glmnet-based methods (multi_enet,
#' joinet, stacking and the univariate elastic net) fit on compressed rows:
#' each training set's n rows are replaced by about p + q rows with the same
#' cross-products, built once per run and shared by the methods. The fits
#' agree with those on the original rows to rounding error. Compression
#' applies where it removes enough rows (p + q + 1 at most
#' \code{getOption('isotwas.glmnet_compress_ratio', 0.75)} times the training
#' set size), and never depends on timings, so a rerun gives identical
#' results. It gains most with an optimised BLAS such as OpenBLAS.
#' \code{options(isotwas.glmnet_compress = FALSE)} turns it off;
#' \code{options(isotwas.glmnet_compress = TRUE)} also turns it on with
#' \code{shared_folds = FALSE} and no \code{folds}, and when the methods are
#' called on their own.
#'
#' @return an isotwas_result. Its \code{comparison} table holds each method's
#'   cross-validated R-squared as the method reported it. Selection scores a
#'   model that selected no SNPs as described in Details, so a method can show
#'   the highest R-squared in \code{comparison} without being selected.
#'
#'
#' @export
compute_isotwas <- function(X,
                            Y,
                            gene_exp = NULL,
                            Y.rep,
                            R = NULL,
                            id,
                            omega_est = 'replicates',
                            r2_scale = c('coherent', 'legacy'),
                            omega_nlambda = 10,
                            method = c('mrce_lasso',
                                       'multi_enet',
                                       'joinet',
                                       'spls',
                                       'sgl',
                                       'mtlasso',
                                       'stacking',
                                       'graph_reg',
                                       'univariate'),
                            similarity_matrix = NULL,
                            predict_nlambda = 20,
                            family = 'gaussian',
                            scale = FALSE,
                            alpha = 0.5,
                            nfolds = 5,
                            verbose = TRUE,
                            par = FALSE,
                            n.cores = NULL,
                            tx_names = NULL,
                            seed = NULL,
                            run_all = TRUE,
                            return_all = FALSE,
                            tol.in = 1e-3,
                            maxit.in = 1e3,
                            joinet_nfolds_int = 10,
                            sgl_nlambda = 20,
                            shared_folds = TRUE,
                            folds = NULL) {

  # Start timing
  start_time <- Sys.time()

  # ============================================================================
  # INPUT VALIDATION
  # ============================================================================

  # Check X
  if (!is.matrix(X) && !is.data.frame(X)) {
    stop("X must be a matrix or data.frame of SNP dosages")
  }
  X <- as.matrix(X)

  # Check Y
  if (!is.matrix(Y) && !is.data.frame(Y)) {
    stop("Y must be a matrix or data.frame of isoform expression")
  }
  Y <- as.matrix(Y)

  # Check dimensions match
  if (nrow(X) != nrow(Y)) {
    stop("Dimension mismatch: X has ", nrow(X), " samples but Y has ", nrow(Y),
         " samples. These must match.")
  }

  # Check Y.rep
  if (missing(Y.rep) || is.null(Y.rep)) {
    stop("Y.rep is required. Provide isoform expression with replicates.")
  }
  Y.rep <- as.matrix(Y.rep)

  # Extract dimensions
  N <- nrow(Y)
  P <- ncol(X)
  G <- ncol(Y)

  # Check/compute R
  if (is.null(R)) {
    R <- nrow(Y.rep) / nrow(Y)
    if (R != floor(R)) {
      stop("Cannot determine R: nrow(Y.rep)=", nrow(Y.rep),
           " is not divisible by nrow(Y)=", nrow(Y),
           ". Please specify R explicitly.")
    }
  } else {
    if (nrow(Y.rep) != R * N) {
      stop("Dimension mismatch: Y.rep has ", nrow(Y.rep), " rows but expected ",
           R * N, " rows (R=", R, " x N=", N, ")")
    }
  }

  # Check id
  if (missing(id) || is.null(id)) {
    stop("id is required. Provide a vector mapping replicates to sample IDs.")
  }

  # Check nfolds
  if (nfolds < 2) {
    stop("nfolds must be at least 2 (got ", nfolds, ")")
  }
  if (nfolds > N) {
    warning("nfolds (", nfolds, ") > N (", N, "). Setting nfolds = N for leave-one-out CV.")
    nfolds <- N
  }

  r2_scale <- match.arg(r2_scale)
  legacy <- identical(r2_scale, 'legacy')
  floor_r2 <- legacy
  empty_r2 <- if (legacy) 0 else -Inf
  # Every method scores R2 as the adjusted R2 of lm(obs ~ pred), as calc.r2()
  # does; legacy scoring keeps the out-of-sample R2 of the five methods that
  # used it
  r2_method <- if (legacy) 'oos' else 'refit'

  # Restore the caller's RNG state on exit
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    caller_seed <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", caller_seed, envir = globalenv()),
            add = TRUE)
  }

  # Set seed
  if (is.null(seed)) {
    seed <- as.integer(Sys.time()) %% 100000L
  }
  set.seed(seed)

  # One set of cross-validation folds for every method that takes them
  cv_folds <- NULL
  if (!is.null(folds)) {
    cv_folds <- .as_folds(folds, N)
  } else if (shared_folds) {
    cv_folds <- .isotwas_folds(N, nfolds, seed)
  }
  if (!is.null(cv_folds)) {
    nfolds <- cv_folds$nfolds
  }

  # The glmnet-based methods share compressed designs built for the same
  # data and rows within this run, and with shared folds they compress by
  # default, since each design then serves several methods (see
  # .glmnet_memo_start)
  .glmnet_memo_start(compress = !is.null(cv_folds))
  on.exit(.glmnet_memo_end(), add = TRUE)

  # Set transcript names
  if (is.null(tx_names)) {
    if (!is.null(colnames(Y))) {
      tx_names <- colnames(Y)
    } else {
      tx_names <- paste0("Transcript_", 1:G)
    }
  }

  # ============================================================================
  # DATA SUMMARY
  # ============================================================================

  if (verbose) {
    message("\n", paste(rep("=", 60), collapse = ""))
    message("                    isoTWAS Model Training")
    message(paste(rep("=", 60), collapse = ""))
    message("\nData summary:")
    message("  Samples (N):           ", N)
    message("  SNPs (P):              ", format(P, big.mark = ","))
    message("  Transcripts (G):       ", G)
    message("  Replicates (R):        ", R)
    message("  CV folds:              ", nfolds)
    message("  Random seed:           ", seed)

    # Warnings for potential issues
    if (N < 50) {
      message("\n[!] Note: Small sample size (N=", N, "). Results may be unstable.")
    }
    if (P > 10 * N) {
      message("[!] Note: High-dimensional setting (P >> N). Regularization is important.")
    }
    if (G > 20) {
      message("[!] Note: Many transcripts (G=", G, "). This may take a while.")
    }
  }

  # ============================================================================
  # HANDLE SINGLE TRANSCRIPT CASE
  # ============================================================================

  if (G == 1) {
    if (verbose) {
      message("\n[!] Single transcript detected. Using univariate methods only.")
    }
    method <- c('univariate')
    run_all <- FALSE

    if (scale) {
      warning("scale = TRUE is ignored for a single transcript: ",
              "no isoform covariance is estimated when G = 1.")
      scale <- FALSE
    }
  }

  # ============================================================================
  # DETERMINE METHODS TO RUN
  # ============================================================================

  all_methods <- isotwas_methods()

  if (run_all) {
    method <- all_methods
  } else {
    # Validate method names
    invalid <- setdiff(method, all_methods)
    if (length(invalid) > 0) {
      stop("Unknown method(s): ", paste(invalid, collapse = ", "),
           "\nAvailable methods: ", paste(all_methods, collapse = ", "))
    }
  }

  # Handle graph_reg without similarity matrix
  if ('graph_reg' %in% method && is.null(similarity_matrix)) {
    method <- setdiff(method, 'graph_reg')
    if (verbose) {
      message("[!] Note: graph_reg removed (no similarity_matrix provided)")
    }
  }

  # Validate similarity matrix dimensions if provided
  if (!is.null(similarity_matrix)) {
    if (!is.matrix(similarity_matrix)) {
      stop("similarity_matrix must be a matrix")
    }
    if (nrow(similarity_matrix) != G || ncol(similarity_matrix) != G) {
      stop("similarity_matrix dimensions (", nrow(similarity_matrix), "x",
           ncol(similarity_matrix), ") must match number of transcripts (", G, "x", G, ")")
    }
  }

  n_methods <- length(method)

  if (verbose) {
    message("\nMethods to run (", n_methods, "):")
    for (m in method) {
      message("  - ", m)
    }
    message("")
  }

  # ============================================================================
  # COMPUTE OMEGA (error covariance)
  # ============================================================================

  if (verbose) {
    message(paste(rep("-", 60), collapse = ""))
    message("Step 1/3: Computing error covariance (Omega)")
    message(paste(rep("-", 60), collapse = ""))
  }

  omega_list <- NULL
  if (G > 1) {
    omega_list <- compute_omega(Y,
                                Y.rep,
                                R,
                                id,
                                method = omega_est,
                                nlambda = omega_nlambda,
                                verbose = verbose)
    if (verbose) {
      message("  [OK] Omega estimation complete\n")
    }
  }

  # ============================================================================
  # INITIALIZE RESULTS
  # ============================================================================

  r2_mat <- data.frame(Transcript = tx_names)
  r2_mat <- cbind(r2_mat, matrix(NA, nrow = G, ncol = length(all_methods)))
  colnames(r2_mat)[-1] <- all_methods

  all_models <- list()
  method_counter <- 0
  failed_methods <- character(0)
  failed_messages <- character(0)

  # ============================================================================
  # RUN PREDICTION METHODS
  # ============================================================================

  if (verbose) {
    message(paste(rep("-", 60), collapse = ""))
    message("Step 2/3: Training prediction models")
    message(paste(rep("-", 60), collapse = ""))
  }

  # Helper function to run a method safely
  run_method <- function(method_name, method_func) {
    method_counter <<- method_counter + 1

    if (verbose) {
      message("\n[", method_counter, "/", n_methods, "] ", method_name)
    }

    method_start <- Sys.time()

    result <- tryCatch({
      method_func()
    }, error = function(e) {
      if (verbose) {
        message("      [ERROR] ", conditionMessage(e))
      }
      failed_methods <<- c(failed_methods, method_name)
      failed_messages <<- c(failed_messages, conditionMessage(e))
      return(NULL)
    })

    if (!is.null(result)) {
      method_time <- round(difftime(Sys.time(), method_start, units = "secs"), 1)
      mean_r2 <- mean(result$summary$r2, na.rm = TRUE)

      if (verbose) {
        message("      [OK] Complete (", method_time, "s) | Mean R2: ", sprintf("%.3f", mean_r2))
      }
    }

    return(result)
  }

  # --- SPLS ---
  if ('spls' %in% method) {
    spls_mod <- run_method("Sparse PLS (spls)", function() {
      multivariate_spls(X = X,
                        Y = Y,
                        nfolds = nfolds,
                        verbose = FALSE,
                        tx_names = tx_names,
                        par = par,
                        n.cores = n.cores,
                        seed = seed,
                        folds = cv_folds)
    })
    if (!is.null(spls_mod)) {
      all_models <- rlist::list.append(all_models, spls_mod)
      r2_mat[, 'spls'] <- spls_mod$summary$r2
    }
  }

  # --- MRCE LASSO ---
  if ('mrce_lasso' %in% method) {
    mrce_mod <- run_method("MRCE Lasso (mrce_lasso)", function() {
      # The folds and cross-products are the same for every Omega
      mrce_cache <- .mrce_fold_cache(X, Y, nfolds, seed, cv_folds)
      run_mrce_omega <- function(omega_idx) {
        compute_mrce(X = X,
                     Y = Y,
                     lambda = NULL,
                     nlambda = predict_nlambda,
                     Omega = omega_list$icov[[omega_idx]],
                     nfolds = nfolds,
                     tol.in = tol.in,
                     maxit.in = maxit.in,
                     verbose = FALSE,
                     seed = seed,
                     parallel = FALSE,
                     n.cores = 1,
                     r2_method = r2_method,
                     floor_r2 = floor_r2,
                     fold_cache = mrce_cache,
                     folds = cv_folds)
      }

      if (par && !is.null(n.cores) && n.cores > 1) {
        mrce_lasso_list <- .mclapply_seeded(
          1:omega_nlambda,
          run_mrce_omega,
          mc.cores = min(n.cores, omega_nlambda),
          seed = seed
        )
      } else {
        mrce_lasso_list <- lapply(1:omega_nlambda, run_mrce_omega)
      }

      get_best(mrce_lasso_list, G = G, empty_r2 = empty_r2)
    })
    if (!is.null(mrce_mod)) {
      all_models <- rlist::list.append(all_models, mrce_mod)
      r2_mat[, 'mrce_lasso'] <- mrce_mod$summary$r2
    }
  }

  # --- MULTIVARIATE ELASTIC NET ---
  if ('multi_enet' %in% method) {
    enet_mod <- run_method("Multivariate Elastic Net (multi_enet)", function() {
      multivariate_elasticnet(X = X,
                              Y = Y,
                              Omega = omega_list$icov[[omega_nlambda]],
                              scale = scale,
                              alpha = alpha,
                              nfolds = nfolds,
                              verbose = FALSE,
                              par = par,
                              n.cores = n.cores,
                              tx_names = tx_names,
                              seed = seed,
                              folds = cv_folds)
    })
    if (!is.null(enet_mod)) {
      all_models <- rlist::list.append(all_models, enet_mod)
      r2_mat[, 'multi_enet'] <- enet_mod$summary$r2
    }
  }

  # --- JOINET ---
  if ('joinet' %in% method) {
    joinet_mod <- run_method("Joinet (joinet)", function() {
      multivariate_joinet(X = X,
                          Y = Y,
                          Omega = omega_list$icov[[omega_nlambda]],
                          scale = scale,
                          alpha = alpha,
                          nfolds = nfolds,
                          verbose = FALSE,
                          par = par,
                          n.cores = n.cores,
                          tx_names = tx_names,
                          seed = seed,
                          nfolds.int = joinet_nfolds_int,
                          folds = cv_folds)
    })
    if (!is.null(joinet_mod)) {
      all_models <- rlist::list.append(all_models, joinet_mod)
      r2_mat[, 'joinet'] <- joinet_mod$summary$r2
    }
  }

  # --- SPARSE GROUP LASSO ---
  if ('sgl' %in% method) {
    sgl_mod <- run_method("Sparse Group Lasso (sgl)", function() {
      multivariate_sgl(X = X,
                       Y = Y,
                       alpha = alpha,
                       nlambda = sgl_nlambda,
                       nfolds = nfolds,
                       verbose = FALSE,
                       seed = seed,
                       par = par,
                       n.cores = n.cores,
                       r2_method = r2_method,
                       floor_r2 = floor_r2,
                       folds = cv_folds)
    })
    if (!is.null(sgl_mod)) {
      all_models <- rlist::list.append(all_models, sgl_mod)
      r2_mat[, 'sgl'] <- sgl_mod$summary$r2
    }
  }

  # --- MULTI-TASK LASSO ---
  if ('mtlasso' %in% method) {
    mtlasso_mod <- run_method("Multi-task Lasso (mtlasso)", function() {
      multivariate_mtlasso(X = X,
                           Y = Y,
                           regularization = "L21",
                           nfolds = nfolds,
                           verbose = FALSE,
                           seed = seed,
                           par = par,
                           n.cores = n.cores,
                           r2_method = r2_method,
                           floor_r2 = floor_r2,
                           folds = cv_folds)
    })
    if (!is.null(mtlasso_mod)) {
      all_models <- rlist::list.append(all_models, mtlasso_mod)
      r2_mat[, 'mtlasso'] <- mtlasso_mod$summary$r2
    }
  }

  # --- STACKING ---
  if ('stacking' %in% method) {
    stacking_mod <- run_method("Stacking (stacking)", function() {
      multivariate_stacking(X = X,
                            Y = Y,
                            Y.rep = Y.rep,
                            R = R,
                            id = id,
                            Omega = omega_list$icov[[omega_nlambda]],
                            nfolds_stack = nfolds,
                            verbose = FALSE,
                            seed = seed,
                            parallel = par,
                            n.cores = n.cores,
                            r2_method = r2_method,
                            floor_r2 = floor_r2,
                            folds = cv_folds)
    })
    if (!is.null(stacking_mod)) {
      all_models <- rlist::list.append(all_models, stacking_mod)
      r2_mat[, 'stacking'] <- stacking_mod$summary$r2
    }
  }

  # --- GRAPH-REGULARIZED ---
  if ('graph_reg' %in% method) {
    graph_mod <- run_method("Graph-regularized (graph_reg)", function() {
      multivariate_graph_reg(X = X,
                             Y = Y,
                             similarity_matrix = similarity_matrix,
                             alpha = alpha,
                             nfolds = nfolds,
                             verbose = FALSE,
                             seed = seed,
                             par = par,
                             n.cores = n.cores,
                             r2_method = r2_method,
                             floor_r2 = floor_r2,
                             folds = cv_folds)
    })
    if (!is.null(graph_mod)) {
      all_models <- rlist::list.append(all_models, graph_mod)
      r2_mat[, 'graph_reg'] <- graph_mod$summary$r2
    }
  }

  # --- UNIVARIATE ---
  if ('univariate' %in% method) {
    uni_mod <- run_method("Univariate ensemble (univariate)", function() {
      uni_enet <- univariate_elasticnet(X = X,
                                        Y = Y,
                                        Omega = omega_list$icov[[omega_nlambda]],
                                        family = family,
                                        scale = scale,
                                        alpha = alpha,
                                        nfolds = nfolds,
                                        verbose = FALSE,
                                        par = par,
                                        n.cores = n.cores,
                                        tx_names = tx_names,
                                        seed = seed,
                                        folds = cv_folds)

      uni_blup <- univariate_blup(X = X,
                                  Y = Y,
                                  Omega = omega_list$icov[[omega_nlambda]],
                                  scale = scale,
                                  nfolds = nfolds,
                                  verbose = FALSE,
                                  par = par,
                                  n.cores = n.cores,
                                  tx_names = tx_names,
                                  seed = seed,
                                  folds = cv_folds)

      uni_susie <- univariate_susie(X = X,
                                    Y = Y,
                                    Omega = omega_list$icov[[omega_nlambda]],
                                    scale = scale,
                                    alpha = alpha,
                                    nfolds = nfolds,
                                    verbose = FALSE,
                                    par = par,
                                    n.cores = n.cores,
                                    tx_names = tx_names,
                                    seed = seed,
                                    folds = cv_folds)

      get_best(list(uni_enet, uni_blup, uni_susie), G = G,
               empty_r2 = empty_r2)
    })
    if (!is.null(uni_mod)) {
      all_models <- rlist::list.append(all_models, uni_mod)
      r2_mat[, 'univariate'] <- uni_mod$summary$r2
    }
  }

  # Check if any models succeeded
  if (length(all_models) == 0) {
    stop("All methods failed. Please check your data and parameters.")
  }

  # Remove all-NA columns from r2_mat
  r2_mat <- r2_mat[, apply(r2_mat, 2, function(x) !all(is.na(x))), drop = FALSE]

  # ============================================================================
  # SELECT BEST MODEL
  # ============================================================================

  if (verbose) {
    message("\n", paste(rep("-", 60), collapse = ""))
    message("Step 3/3: Selecting best model")
    message(paste(rep("-", 60), collapse = ""))
  }

  isotwas_mod <- get_best(all_models, G = G, empty_r2 = empty_r2)
  colnames(Y) <- tx_names

  # ============================================================================
  # GENE EXPRESSION PREDICTION (optional)
  # ============================================================================

  tx2gene_coef <- NULL

  if (!is.null(gene_exp)) {
    if (verbose) {
      message("\nFitting transcript-to-gene aggregation model...")
    }

    pred_mat <- get_prediction_matrix(isotwas_mod)

    set.seed(seed)
    test.folds <- caret::createFolds(1:nrow(Y),
                                     k = nfolds,
                                     returnTrain = FALSE,
                                     list = FALSE)
    if (!is.null(cv_folds)) {
      test.folds <- cv_folds$foldid
    }

    glmnet_pred <- glmnet::cv.glmnet(x = pred_mat,
                                     y = gene_exp,
                                     foldid = test.folds,
                                     keep = TRUE,
                                     relax = TRUE,
                                     gamma = c(0, .25, .5, .75, 1),
                                     trace.it = FALSE)

    gamma_which <- which.max(sapply(lapply(glmnet_pred$fit.preval,
                                           function(y) {
                                             apply(y, 2, function(x) {
                                               summary(stats::lm(gene_exp ~ x))$adj.r.sq
                                             })
                                           }), max))

    which_r2 <- which.max(apply(glmnet_pred$fit.preval[[gamma_which]],
                                2,
                                function(x) {
                                  summary(stats::lm(gene_exp ~ x))$adj.r.sq
                                }))

    tot_mod <- glmnet::glmnet(x = pred_mat,
                              y = gene_exp,
                              alpha = c(0, .25, .5, .75, 1)[gamma_which],
                              lambda = glmnet_pred$lambda[which_r2])

    glmnet_r2 <- max(sapply(lapply(glmnet_pred$fit.preval,
                                   function(y) {
                                     apply(y, 2, function(x) {
                                       summary(stats::lm(gene_exp ~ x))$adj.r.sq
                                     })
                                   }), max))

    if (all(as.numeric(coef(tot_mod))[-1] == 0)) {
      ccc <- stats::coef(stats::lm(gene_exp ~ pred_mat))[-1]
    } else {
      ccc <- as.numeric(stats::coef(tot_mod))[-1]
    }

    tx2gene_coef <- data.frame(Feature = tx_names,
                               Weight_tx2gene = ccc,
                               R2 = glmnet_r2)

    if (verbose) {
      message("  [OK] Transcript-to-gene R2: ", sprintf("%.3f", glmnet_r2))
    }
  }

  # ============================================================================
  # BUILD RESULT
  # ============================================================================

  result <- create_isotwas_result(
    best_models = isotwas_mod,
    method_results = if (return_all) all_models else NULL,
    comparison_df = r2_mat,
    tx2gene = tx2gene_coef,
    failed_methods = data.frame(method = failed_methods,
                                message = failed_messages,
                                stringsAsFactors = FALSE),
    seed = seed
  )

  # ============================================================================
  # FINAL SUMMARY
  # ============================================================================

  total_time <- round(difftime(Sys.time(), start_time, units = "secs"), 1)

  if (verbose) {
    message("\n", paste(rep("=", 60), collapse = ""))
    message("                         RESULTS")
    message(paste(rep("=", 60), collapse = ""))

    # Method comparison
    message("\nMethod comparison (mean R2 across transcripts):")
    method_means <- sapply(colnames(r2_mat)[-1], function(m) {
      mean(r2_mat[[m]], na.rm = TRUE)
    })
    method_means <- sort(method_means, decreasing = TRUE)

    for (i in seq_along(method_means)) {
      marker <- if (i == 1) " <-- BEST" else ""
      message("  ", sprintf("%-15s", names(method_means)[i]),
              sprintf("%.4f", method_means[i]), marker)
    }

    # Best model info
    best_method <- isotwas_mod$method
    best_mean_r2 <- mean(isotwas_mod$summary$r2, na.rm = TRUE)

    message("\nBest model: ", best_method)
    message("  Mean R2:  ", sprintf("%.4f", best_mean_r2))
    message("  R2 range: [", sprintf("%.4f", min(isotwas_mod$summary$r2, na.rm = TRUE)),
            ", ", sprintf("%.4f", max(isotwas_mod$summary$r2, na.rm = TRUE)), "]")

    # Per-transcript summary
    if (G <= 10) {
      message("\nPer-transcript R2 (best model):")
      for (i in 1:G) {
        message("  ", sprintf("%-20s", tx_names[i]),
                sprintf("%.4f", isotwas_mod$summary$r2[i]))
      }
    }

    # Warnings
    if (length(failed_methods) > 0) {
      message("\n[!] Warning: ", length(failed_methods), " method(s) failed: ",
              paste(failed_methods, collapse = ", "))
    }

    if (best_mean_r2 < 0.01) {
      message("\n[!] Warning: Very low R2. The model may not be predictive.")
    }

    message("\nTotal time: ", total_time, " seconds")
    message(paste(rep("=", 60), collapse = ""), "\n")
  }

  return(result)
}

