#' Compute isoTWAS model for a set of isoforms/transcripts
#'
#' The function runs MVR models for a set of transcripts and outputs
#' the best model
#'
#' @param X matrix, design matrix of SNP dosages
#' @param Y matrix, matrix of G isoform expression across columns
#' @param Y.rep matrix, matrix of G isoform expression with replicates
#' @param R int, number of replicates
#' @param id vector, vector of sample ids showing rep to id
#' @param omega_est character, 'replicates' or 'mean' to use Y.rep or Y
#' @param omega_nlambda int, number of omegas to generate
#' @param method character, vector of methods to use
#' @param predict_nlambda int, number of lambdas in MRCE
#' @param family character, glmnet family
#' @param scale logical, T/F to scale Y by Omega
#' @param alpha numeric, elastic net mixing parameter
#' @param nfolds int, number of CV folds
#' @param verbose logical
#' @param par logical, uses mclapply to parallelize model fit
#' @param n.cores int, number of parallel cores
#' @param tx_names vector, character vector of tx names - order of columns of Y
#' @param seed int, random seed
#' @param return_all logical, return R2 for all models?
#' @param tol.in numeric, tolerance for objective difference
#' @param maxit.in int, maximum number of iteractions
#'
#' @return optimal isoTWAS model
#'
#'
#' @export
compute_isotwas <- function(X,
                            Y,
                            gene_exp = NULL,
                            Y.rep,
                            R,
                            id,
                            omega_est = 'replicates',
                            omega_nlambda = 10,
                            method = c('mrce_lasso',
                                       'curds_whey',
                                       'multi_enet',
                                       'joinet',
                                       'mvsusie',
                                       'mrmash',
                                       'finemap',
                                       'univariate'),
                            predict_nlambda = 50,
                            family = 'gaussian',
                            scale = F,
                            alpha = 0.5,
                            nfolds = 5,
                            verbose = F,
                            par = F,
                            n.cores = NULL,
                            tx_names = NULL,
                            seed = NULL,
                            run_all = T,
                            return_all = F,
                            tol.in = 1e-3,
                            maxit.in = 1e3,
                            coverage = .9){

  ### CHECKS
  if (nrow(X) != nrow(Y)){
    stop('No. of rows of X =/= no. of rows of Y.')
  }

  N = nrow(Y)
  P = ncol(X)
  G = ncol(Y)
  pred_mat = matrix(nrow = N,
                    ncol = G)

  if (is.null(R)){
    R = nrow(Y.rep)/nrow(Y)
  } else {
    if (nrow(Y.rep) != R*N){
      stop('No. of rows of Y.full =/= R*N')
    }
  }

  ### COMPUTE OMEGA
  if (verbose){print('Computing Omega')}
  if (ncol(Y) > 1 & ncol(Y.rep) > 1){
    omega_list = compute_omega(Y,
                               Y.rep,
                               R,
                               id,
                               method = omega_est,
                               nlambda = omega_nlambda,
                               verbose = verbose)}
  if (ncol(Y) == 1){

    method = c('univariate')
    Y = as.matrix

  }

  if (run_all){
    print('run_all is set to T and all methods will be run.')
    method = c('mrce_lasso',
               'curds_whey',
               'multi_enet',
               'joinet',
               'mvsusie',
               'mrmash',
               'finemap',
               'univariate')
  }

  if (is.null(tx_names)){
    tx_names = colnames(Y)
  }

  r2_mat = data.frame(Transcript = tx_names)
  r2_mat = cbind(r2_mat,
                 matrix(nrow = nrow(r2_mat),
                        ncol = 8))
  colnames(r2_mat)[-1] =  c('mrce_lasso',
                            'curds_whey',
                            'multi_enet',
                            'joinet',
                            'mvsusie',
                            'mrmash',
                            'finemap',
                            'univariate')

  if (is.null(seed)){
    seed = sample(1:100000,1)
  }

  all_models = list()

  if ('mrce_lasso' %in% method){
    if (verbose){print('Running mrce_lasso')}
      mrce_lasso = list()
      for (i in 1:omega_nlambda){
        mrce_lasso = rlist::list.append(mrce_lasso,
                                        compute_mrce(X = X,
                                                     Y = Y,
                                                     lambda = NULL,
                                                     nlambda =
                                                       predict_nlambda,
                                                     Omega =
                                                       omega_list$icov[[i]],
                                                     nfolds = nfolds,
                                                     tol.in = tol.in,
                                                     maxit.in = maxit.in,
                                                     verbose = verbose,
                                                     seed = seed))
      }
    all_models = rlist::list.append(all_models,get_best(mrce_lasso,
                                                        G = G))
    r2_mat[,'mrce_lasso'] = sapply(mrce_lasso[[1]],function(x) x$R2)
  }

  ## CURDS AND WHEY
  if ('curds_whey' %in% method){
    if (verbose){print('Running curds_whey')}
    best_curds_whey = compute_curds_whey(X,
                                         Y,
                                         family = family,
                                         alpha = alpha,
                                         nfolds = nfolds,
                                         verbose = verbose,
                                         par = F,
                                         n.cores = NULL,
                                         tx_names = tx_names,
                                         seed)
    all_models = rlist::list.append(all_models,best_curds_whey)
    r2_mat[,'curds_whey'] = sapply(best_curds_whey[[1]],
                                   function(x) x$R2)
  }

  ## MULTIVARATE ELASTIC NET
  if ('multi_enet' %in% method){
    if (verbose){print('Running multi_enet')}
    best_multi_enet =
      multivariate_elasticnet(X = X,
                              Y = Y,
                              Omega =
                                omega_list$icov[[omega_nlambda]],
                              scale = scale,
                              alpha = alpha,
                              nfolds = nfolds,
                              verbose = verbose,
                              par = par,
                              n.cores = n.cores,
                              tx_names = tx_names,
                              seed = seed)
    all_models = rlist::list.append(all_models,best_multi_enet)
    r2_mat[,'multi_enet'] = sapply(best_multi_enet,
                                   function(x) x$R2)
  }

  ## JOINET
  if ('joinet' %in% method){
    if (verbose){print('Running joinet')}
    best_joinet =
      multivariate_joinet(X = X,
                          Y = Y,
                          Omega =
                            omega_list$icov[[omega_nlambda]],
                          scale = scale,
                          alpha = alpha,
                          nfolds = nfolds,
                          verbose = verbose,
                          par = par,
                          n.cores = n.cores,
                          tx_names = tx_names,
                          seed = seed)
    all_models = rlist::list.append(all_models,best_joinet)
    r2_mat[,'joinet'] = sapply(best_joinet,
                                   function(x) x$R2)
  }

  # MVSUSIE
  if ('mvsusie' %in% method){
    if (verbose){print('Running mvsusie')}
    mmbr_mod = multivariate_mmbr(X = X,
                                 Y = Y,
                                 nfolds = nfolds,
                                 verbose = verbose,
                                 tx_names = tx_names,
                                 seed = seed,
                                 par = par,
                                 n.cores = n.cores)
    all_models = rlist::list.append(all_models,mmbr_mod)
    r2_mat[,'mvsusie'] = sapply(mmbr_mod,
                               function(x) x$R2)
  }

  # MR MASH
  if ('mrmash' %in% method){
    if (verbose){print('Running mrmash')}
    mrmash_mod = multivariate_mrmash(X = X,
                        Y = Y,
                        nfolds = nfolds,
                        verbose = verbose,
                        tx_names = tx_names,
                        par = par,
                        n.cores = n.cores,
                        seed = seed)
    all_models = rlist::list.append(all_models,mrmash_mod)
    r2_mat[,'mrmash'] = sapply(mrmash_mod,
                               function(x) x$R2)
  }

  if ('finemap' %in% method){
    if (verbose){print('Running finemap')}
    best_finemap = compute_finemap_regress(X = X,
                                           Y = Y,
                                           Y.rep = Y.rep,
                                           R = R,
                                           id = id,
                                           nfolds = nfolds,
                                           verbose = verbose,
                                           tx_names = tx_names,
                                           coverage = coverage,
                                           seed = seed)
    all_models = rlist::list.append(all_models,best_finemap)
    r2_mat[,'finemap'] = sapply(best_finemap,
                               function(x) x$R2)
  }

  ### UNIVARIATE FUSION
  if ('univariate' %in% method){
    if (verbose){print('Running univariate')}
    uni_enet = univariate_elasticnet(X = X,
                                     Y = Y,
                                     Omega = omega_list[[omega_nlambda]],
                                     family = family,
                                     scale = scale,
                                     alpha = alpha,
                                     nfolds = nfolds,
                                     verbose = verbose,
                                     par = F,
                                     n.cores = 1,
                                     tx_names = tx_names,
                                     seed = seed)

    uni_blup = univariate_blup(X = X,
                               Y = Y,
                               Omega = omega_list[[omega_nlambda]],
                               scale = scale,
                               alpha = alpha,
                               nfolds = nfolds,
                               verbose = verbose,
                               par = F,
                               n.cores = 1,
                               tx_names = tx_names,
                               seed = seed)

    uni_susie = univariate_susie(X = X,
                                 Y = Y,
                                 Omega = omega_list[[omega_nlambda]],
                                 scale = scale,
                                 alpha = alpha,
                                 nfolds = nfolds,
                                 verbose = verbose,
                                 par = F,
                                 n.cores = 1,
                                 tx_names = tx_names,
                                 seed = seed)

    univariate = list(uni_enet,
                      uni_blup,
                      uni_susie)
    best_uni = get_best(univariate,
                        G = G)
    all_models = rlist::list.append(all_models,best_uni)
    r2_mat[,'univariate'] = sapply(best_uni,
                                function(x) x$R2)

  }
  r2_mat = r2_mat[, apply(r2_mat,2,function(x) !all(is.na(x)))]

  isotwas_mod = get_best(all_models,
                         G = G)
  colnames(Y) = tx_names



  if (!is.null(gene_exp)){

    for (i in 1:length(isotwas_mod)){
      pred_mat[,i] = isotwas_mod[[i]]$Pred
    }

    set.seed(seed)
    test.folds = caret::createFolds(1:nrow(Y),
                                     k = nfolds,
                                     returnTrain = F,
                                     list = F)
    glmnet_pred = glmnet::cv.glmnet(x = pred_mat,
                                    y = gene_exp,
                                    foldid = test.folds,
                                    keep = T,
                                    relax = T,
                                    gamma = c(0,.25,.5,.75,1),
                                    trace.it = F)

    gamma_which = which.max(sapply(lapply(glmnet_pred$fit.preval,
           function(y){
             apply(y,2,function(x){
               summary(lm(gene_exp ~ x))$adj.r.sq
             })
           }),max))

    which_r2 = which.max(apply(glmnet_pred$fit.preval[[gamma_which]],
                    2,
                    function(x){summary(lm(gene_exp ~ x))$adj.r.sq}
                    ))

    tot_mod = glmnet(x = pred_mat,
                     y = gene_exp,
                     alpha = c(0,.25,.5,.75,1)[gamma_which],
                     lambda = glmnet_pred$lambda[which_r2])
    glmnet_r2 = max(sapply(lapply(glmnet_pred$fit.preval,
                                  function(y){
                                    apply(y,2,
                                          function(x){
                                            summary(lm(gene_exp ~ x))$adj.r.sq
                                          })
                                  }),max))

    if (all(as.numeric(coef(tot_mod))[-1] == 0)){
      ccc = coef(lm(gene_exp ~ pred_mat))[-1]
    } else {
      ccc = as.numeric(coef(tot_mod))[-1]
    }
    tx2gene_coef = data.frame(Feature = tx_names,
                              Weight_tx2gene =
                                ccc,
                              R2 = glmnet_r2)


  } else {
    tx2gene_coef = 'Gene expression vector not supplied'
  }

  if (return_all){

    r2.df = r2_mat
    isotwas_mod = list(Model = get_best(all_models,
                                        G = G),
                       R2 = r2.df)
  }

  return(list(isotwas_mod = isotwas_mod,
              tx2gene_coef = tx2gene_coef))


}
