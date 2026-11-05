# isotwas 2.0.0.9000

## Faster model training, same estimates

Seven fitting routines are reimplemented to do the same computation with
less work. Each reproduces what it replaces: the same folds and random
numbers, tuning grids, update order and stopping rules. Fitted weights,
cross-validated predictions and R-squared agree to rounding error except
where noted. `tests/testthat/test-fast-*.R` check each against the
original. `benchmarks/` (not installed) has the timing script, and the
speed-ups below were measured with it (OpenBLAS, one core, simulated
genes).

* The glmnet-based methods (`multivariate_elasticnet()`,
  `univariate_elasticnet()`, `multivariate_joinet()`, `multivariate_stacking()`)
  can fit each elastic-net path on compressed rows. A gaussian or mgaussian
  path depends on the data only through the centred X'X, X'Y and Y'Y, so
  `R/glmnet_compress.R` replaces the n rows by rank(X) + 1 rows per response
  (rank(X) + q for mgaussian) with the same cross-products, from a pivoted
  Cholesky factor of X'X, and calls `glmnet()` on those. Each training fold's
  cross-products are the full ones minus its held-out rows. `.cv_glmnet()`
  reproduces `cv.glmnet(keep = TRUE)` around these fits with the same folds,
  lambda interpolation and choice of `lambda.min`. With `verbose = TRUE`,
  `multivariate_elasticnet()` and `multivariate_joinet()` no longer print
  glmnet's progress bars, and `par = TRUE` runs their folds with `mclapply()`
  rather than a doParallel backend, as `univariate_elasticnet()` already did
  (so serially on Windows). joinet is reimplemented in
  R (`R/joinet_fast.R`) the same way, with joinet 1.0.0's random draws, and no
  longer fits an unused mgaussian model in every external fold; joinet moves
  to Suggests. Results agree with the originals to 1e-11 to 1e-14, and are
  identical where the original rows are used.

  Compression is off unless `options(isotwas.glmnet_compress = TRUE)`. It
  then applies wherever it removes enough rows (p + q + 1 at most
  `getOption("isotwas.glmnet_compress_ratio", 0.75)` times the training fold
  size), never depending on timings, so a rerun gives identical results. It
  saves glmnet time roughly in proportion to the rows removed, while building
  the rows is BLAS work, so whether it is faster depends on the BLAS. With
  OpenBLAS on the cluster, at n = 2011, the four methods ran 1.2-2.3x as fast
  at p = 550-850 and 0.9-1.5x at p = 1150. With R's reference BLAS, building
  the rows can cost more than it saves (at p = 850: 0.5-1.3x, slowest for
  `univariate_elasticnet()`, whose fold designs each serve one transcript).

* `univariate_blup()` no longer calls `rrBLUP::mixed.solve()` once per
  transcript. Its kernel is a function of X'X, so both spectra REML needs come
  from one eigendecomposition of the smaller Gram matrix (X'X or XX') and one
  of a diagonal-minus-rank-one matrix per training set, shared by every
  transcript. `mixed.solve()` did two SVDs, a complete QR and an n x n inverse
  per transcript. The REML objective, bounds and `optimize()` call are
  unchanged. Weights agree with `mixed.solve()` to about 1e-5 relative:
  `optimize()` locates the REML optimum only to about 1e-4 relative, and
  `mixed.solve()` itself moves by up to 2e-6 when its input moves by one ulp.
  Missing values fall back to `mixed.solve()`. 4.5-34x faster, more at larger n.

* `univariate_susie()` no longer calls `susieR::susie()` on a sparse copy of
  the dosages. `.susie_fit()` runs the same IBSS iterations as susieR 0.14.2
  (the CRAN release) for this fixed-prior-variance configuration, whatever
  susieR version is installed (see Changes to fitted models). It forms X'r
  from X'y and the training set's standardised X'X (or from cached fitted
  values when p > n). The training sets' X'X come from one full crossprod
  minus the held-out rows. Agrees with susieR 0.14.2's `susie()` to 1e-15,
  iteration count included. 2.8-7.5x faster.

* `multivariate_spls()` no longer calls `spls::cv.spls()` and `spls::spls()`.
  `.cv_spls()` and `.spls_fit()` reproduce them for the arguments
  multivariate_spls uses, calling `pls::simpls.fit()` directly instead of
  through a model frame, centring each training set once for all nine `eta`,
  and forming the direction vectors' Z Z' v products as Z (Z' v). Results are
  bit-identical in every benchmark. 2.4-10x faster, more at larger p. `pls`
  is now in Imports.

* `multivariate_sgl()`, `multivariate_mtlasso()` (L21 backend) and
  `multivariate_graph_reg()` share one compiled block coordinate-descent solver
  (`src/rowcd.c`). It makes the same row updates in the same order as the R
  solvers, which are kept as `.fit_sgl_r()`, `.fit_mtl_l21_r()` and
  `.fit_graph_reg_r()`. The gradient comes from X'X, updated by one column when a
  row changes, rather than from three n x q residual updates per SNP.
  graph_reg's per-SNP q x q `solve()` becomes a product with the Laplacian's
  eigenvectors. Each training set's X'X is the full X'X minus its held-out
  rows. Above `getOption("isotwas.gram_max_p", 5000)` SNPs the solver works
  from residuals instead. Agrees with the R solvers to 1e-14. 7-67x faster.

* MRCE's coordinate descent (`rblasso()`) uses a new kernel,
  `src/mrce-fast.c`. `blasso()` updates all p*q soft-threshold values whenever
  one coefficient changes. The new kernel keeps S B instead, which one change
  touches in one column, and forms each soft-threshold value when it is
  visited. Agrees with `blasso()` to 1e-14 with the same iteration counts,
  warm starts included. `options(isotwas.mrce_reference = TRUE)` restores
  the old kernel. `compute_mrce()` gains `fold_cache`, which `compute_isotwas()`
  uses to compute the folds' variance filters and cross-products once for all
  `omega_nlambda` Omegas rather than once per Omega. It also computes the
  cor.test p-value only at the selected lambda, the only one used. 3.4-7.2x
  faster for one Omega, 4.5-7.1x for `compute_isotwas()`'s 10 Omegas.

* `univariate_elasticnet()`, `univariate_blup()` and `univariate_susie()` no
  longer run a `cor.test()` they discard when scoring R-squared.

## Changes to fitted models

* `compute_isotwas()` gives every method the same cross-validation folds
  (`shared_folds = TRUE`, the default), so the cross-validated R-squared that
  `get_best_isotwas()` compares across methods come from the same splits. The
  folds are those drawn by `set.seed(seed)` then
  `caret::createFolds(1:N, k = nfolds, returnTrain = TRUE)`, which
  `univariate_blup()`, `univariate_susie()`, `compute_mrce()`,
  `multivariate_sgl()`, `multivariate_graph_reg()` and `multivariate_spls()`
  already drew for themselves, so their results do not change. The others
  did not: `multivariate_elasticnet()` drew a different fold set for each of
  its three alphas, `univariate_elasticnet()` one per transcript,
  `multivariate_joinet()` its own external folds, `multivariate_stacking()`
  folds drawn after its tuning had moved the random stream, and
  `multivariate_mtlasso()` a second set for scoring R-squared. With the same
  seed their R-squared, and for all but mtlasso their weights, now differ
  from earlier versions. mtlasso's R-squared now comes from the folds that
  chose its lambda, as every other method's does, so it can be slightly
  higher. `shared_folds = FALSE` restores each method's own draw, and
  `folds` takes the user's own folds (training-row indices per fold, or one
  fold label per sample). Each fitting function also gains `folds`; its
  default, `NULL`, keeps the function's own draw.

* `univariate_elasticnet()` and the final fit of `univariate_susie()` now fit
  an intercept. Both used `intercept = FALSE` on expression that still carries
  each feature's mean, so the weights largely reproduced the expression level
  rather than cis-genetic effects: dense, almost all one sign, with a sum that
  tracked the feature mean. Weights are now unchanged by adding a constant to
  expression. SuSiE's cross-validation folds already fitted an intercept, so
  its reported R-squared did not describe the weights it returned.

* A cross-validated prediction that is constant or not positively correlated
  with expression now scores an R-squared of 0, and its model is emptied, as
  for a model that selected no SNPs. `calc.r2()` and the univariate methods
  score by the adjusted R-squared of `lm(obs ~ pred)`, which refits the slope,
  so an anti-correlated prediction used to score a positive R-squared. This
  happened mostly where cross-validation chose the null end of the lambda
  path: multi_enet at `alpha = 0` then returns weights of about 1e-38, which
  were kept and could win selection. Weights all below 1e-8 in absolute value
  are now also treated as empty. multi_enet's choice of `alpha` uses the same
  score, so it changes too.

* A warm-started MRCE lasso fit now converges to the lasso solution.
  `rblasso()` started its soft-threshold values at M - S B0 Omega, which
  still holds each coefficient's own term S[r,r] Omega[c,c] B0[r,c].
  `blasso()` never updates an entry's soft value when that entry changes, so
  the offset stayed for the whole fit: restarting at the exact solution moved
  it, for example from 0.92 to 0 for the largest coefficient. `rblasso()` now
  adds the own term back. Cold starts are unchanged, so the final model at a
  given lambda is the same. The cross-validation fits along the lambda path
  are warm-started, so the chosen lambda and the reported cross-validated
  R-squared can change, in either direction and sometimes by a lot. On 16
  simulated genes (n = 300, p = 150, q = 4) the lambda changed for 4 of them
  and the mean R-squared went from 0.044 to 0.040. On the simulated gene of
  `benchmarks/` (n = 300, p = 400, q = 8) it rose for all 8 isoforms,
  by up to 0.31, with the same weights. Before the fix, the reported R-squared came from fits that were not
  the estimator whose weights were returned.

* `univariate_susie()` no longer depends on the installed susieR version. It
  used `susieR::susie()`'s defaults, and susieR 0.16.x (GitHub) runs IBSS
  differently from 0.14.2 (CRAN): a tolerance of 1e-4 rather than 1e-3, an
  ELBO decrease not counted as converged, prior log odds rather than log
  posterior odds 0 for a SNP with no variance (common within
  cross-validation folds), and trimming of null effects. It now always
  follows 0.14.2, so fits are unchanged for users of the CRAN release and
  change slightly for users of 0.16.x. `finemap_regress()` still calls
  `susieR::susie()` and so still follows the installed version.

* `compute_isotwas()` now defaults to `r2_scale = 'coherent'` (see New
  features). Clamped methods no longer win transcripts that no method
  predicts, and a transcript for which every model is empty is
  labelled `NA` instead of being assigned an empty model. The fitting
  functions' `floor_r2` now defaults to `FALSE` and `get_best_isotwas()`'s
  `empty_r2` to `-Inf`, to match. `r2_scale = 'legacy'` reproduces earlier
  selection.

## Bug fixes

* `compute_omega()` now calls `match.arg()` on `method`. Its documented default
  is `c('replicates','mean')`, which reached `if (method == 'mean')` as a
  length-2 condition — an error since R 4.2, so the exported function could not
  be called at its own default. An unrecognised `method` is now rejected by
  name instead of falling through both branches and returning `NULL`.

* `compute_omega()` releases its output sink with `on.exit()`. An error between
  `sink(tempfile())` and the closing `sink()` previously left the caller's
  console redirected for the rest of the session.

* `univariate_elasticnet(par = TRUE)` no longer fails with "object 'x' not
  found". The parallel branch built its list of responses from an undefined
  lowercase `x`. Because `compute_isotwas()` runs each method inside a
  `tryCatch`, this removed the whole `univariate` method from the comparison
  silently whenever `par = TRUE` and `verbose = FALSE` — and for a
  single-transcript gene, where the method set is forced to `univariate`, it
  meant "All methods failed".

* `compute_isotwas()` passes `omega_list$icov[[omega_nlambda]]` to the
  univariate methods, matching the multivariate ones. They previously received
  `omega_list[[omega_nlambda]]`, which indexes a component of the `huge` object
  itself rather than a precision matrix. This was inert at the default
  `scale = FALSE` and silently wrong at `scale = TRUE`.

* `compute_isotwas()` prunes all-NA columns from the comparison matrix with
  `drop = FALSE`, so a run in which only one method survives still returns a
  matrix rather than a vector.

* `compute_isotwas()` ignores `scale = TRUE` with a warning when `G == 1`. No
  isoform covariance is estimated for a single transcript, so `Omega` is `NULL`
  and scaling reached `eigen(NULL)`.

## Breaking changes

* Four fitting functions now label their output with the name the method
  appears under in the comparison table, so `model$method` can be matched
  against `colnames(result$comparison)` directly:

  | was | now |
  |---|---|
  | `multivariate_enet` | `multi_enet` |
  | `sparse_group_lasso` | `sgl` |
  | `mtlasso_l21` | `mtlasso` |
  | `graph_regularized` | `graph_reg` |

  The RMTL backend of `multivariate_mtlasso()` built its label as
  `paste0("mtlasso_", tolower(regularization))`, so it also emitted
  `mtlasso_trace` and `mtlasso_lasso`; all three are now `mtlasso`. The variant
  is recorded in `result$regularization` and on each transcript model as
  `$regularization`, so it survives into the selected model.

  Code that matches on the old strings needs updating. The univariate leaf
  labels (`univariate_enet`, `univariate_blup`, `univariate_susie`) are
  unchanged: they name members of the `univariate` ensemble rather than
  comparison columns.

## New features

* `compute_omega()` gains `method = 'residuals'`, which runs the graphical
  lasso on the within-sample deviations of `Y.rep` rather than on the fitted
  values of a fixed-effect fit on `id`. The existing `'replicates'` computes
  `H %*% Y.rep`, which replaces every replicate with its sample's mean, so the
  within-sample quantification noise is removed rather than extracted. At
  `R = 1` the design is square, the hat matrix is the identity, and
  `'replicates'` is identical to `'mean'`.

  `'replicates'` remains the default, so existing results reproduce.
  `'residuals'` requires at least two replicates per sample and errors
  otherwise, since the within-sample covariance is unidentifiable at `R = 1`.

* `compute_isotwas()` gains `r2_scale`. Methods are compared by a `which.max`
  over a transcripts-by-methods matrix of cross-validated R-squared, but two
  scales are in use: `mrce_lasso`, `sgl`, `mtlasso`, `graph_reg` and
  `stacking` compute a direct out-of-sample R-squared and clamp it at 0, while
  every other method is scored by `calc.r2()`, the adjusted R-squared of
  `lm(obs ~ pred)`, which passes negatives through. A model that selected no
  SNPs is scored 0 as well.

  On a transcript no method can predict, a clamped method's 0 therefore beats
  every rival's negative score, and an empty model beats a real one.
  `r2_scale = 'coherent'` removes the clamp and scores an empty model `-Inf`,
  leaving a transcript no model could fit labelled `NA` rather than assigned
  an empty winner. `'coherent'` is the default; `'legacy'` reproduces earlier
  selection.

  The underlying `r2_method` and `floor_r2` arguments are also available on
  the individual fitting functions.

## Bug fixes (selection labels)

* `get_best_isotwas()` keeps a nested ensemble's leaf label. `compute_isotwas()`
  builds its `univariate` entry by calling `get_best()` a second time, and that
  entry's own `$method` is the literal `"best"`; stamping it over the inner
  label meant every transcript won by the ensemble reported `"best"` instead of
  `univariate_enet`, `univariate_blup` or `univariate_susie`. Which model wins
  is unchanged -- only what it is called.

## Reproducibility

* `univariate_elasticnet(par = TRUE)` is now reproducible, and independent of
  `n.cores`. Once its parallel branch runs (see the fix under Bug fixes), each
  forked `mclapply()` child draws its CV folds from a different RNG state, so
  repeated runs at the same seed gave different answers. Setting
  `RNGkind("L'Ecuyer-CMRG")` would fix the repeatability but not the dependence
  on core count -- with the default `mc.preschedule = TRUE` the stream advances
  per child rather than per job. Each job is now given its own seed, derived
  from the caller's `seed`, which removes both problems.

  The other parallel paths draw their folds before forking, so they were
  already reproducible: `multivariate_mtlasso()`, `multivariate_stacking()` and
  `multivariate_graph_reg()` give the same results as before, in parallel and
  serially. The per-job seeding is applied at all nine `mclapply()` sites so
  that any random draw inside a worker is covered.

  Serial results are unchanged. For `univariate_elasticnet()`, a parallel run
  still differs from a serial one at the same seed: the serial path draws its
  folds from one stream in order, the parallel path gives each job an
  independent seed. Both are reproducible; they are not the same partition.

* `compute_isotwas()` no longer draws a defaulted `seed` with `sample()`, which
  both consumed the caller's RNG stream and made the drawn seed depend on it.
  The function now restores the caller's RNG state on exit, and the seed the
  run actually used is recorded in the result as `$seed` -- without which a run
  left to the default could not be reproduced at all.

## Performance and tuning

* The MRCE coordinate-descent kernel in `src/mrce-solveb.c` now calls
  `R_CheckUserInterrupt()` once per outer sweep. It is a `.C()` call, so
  without it neither Ctrl-C nor `setTimeLimit()` could stop it: one sweep is
  O((p*q)^2) and `maxit` defaults to 1000, so a single gene could hold a core
  for tens of minutes with no way to bound it from R. Checked per sweep rather
  than in the inner loops, where the call overhead would be material; safe
  because every buffer in the kernel is caller-owned.

* `multivariate_joinet()` exposes `nfolds.int`, previously hardcoded to 10 in
  both its branches, and `compute_isotwas()` plumbs it through as
  `joinet_nfolds_int`. joinet is the most expensive method at moderate
  transcript counts and most of that time is this inner CV.

* `compute_isotwas()` plumbs `sgl_nlambda` through to `multivariate_sgl()`,
  which could not be tuned without editing the package.

* `similarity_from_gtf()` accepts an already-imported `GRanges` in place of a
  path, and restricts its own import to `feature.type = "exon"`. It imported
  the entire annotation on every call, so building matrices for many genes
  re-read the whole file each time. It also now splits exons by transcript
  once rather than rebuilding the same subsets inside the O(K^2) pair loop.

  Defaults are unchanged in all four cases.

## Documentation

* `create_similarity_from_exons()` notes that it compares exon IDENTIFIERS, so
  two isoforms score as dissimilar whenever they share sequence through exons
  with different IDs -- which is what happens at an alternative 5' or 3' end.
  `similarity_from_gtf()` is the base-pair measure, and is now cross-referenced.

### Also

* `compute_omega()` no longer forms the `(R*N) x (R*N)` hat matrix. Because the
  design is one-way on sample id, the projection and its residual are exactly
  "replace by, or subtract, each sample's mean", for balanced and unbalanced
  designs alike. This is O(nrow * ncol) rather than O(nrow^2): at N = 1000 and
  R = 50 the explicit matrix is 18.6 GB. Results are unchanged.

## Improvements

* `compute_isotwas()` records failed methods in its result. `isotwas_result`
  gains a `failed_methods` data frame of the method name and the message it
  failed with, and `print()` reports it. Failures were previously visible only
  as a `verbose` console message.

* `isotwas_methods()` is the single source of truth for the method names, used
  both to build the comparison columns and to check each model's label.
