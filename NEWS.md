# isotwas 2.0.0.9000

## Changes to fitted models

* `univariate_elasticnet()` and the final fit of `univariate_susie()` now fit
  an intercept. Both used `intercept = FALSE` on expression that still carries
  each feature's mean, so the weights largely reproduced the expression level
  rather than cis-genetic effects: dense, almost all one sign, with a sum that
  tracked the feature mean. Weights are now unchanged by adding a constant to
  expression. SuSiE's cross-validation folds already fitted an intercept, so
  its reported R-squared did not describe the weights it returned.

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
  an empty winner. `'legacy'` is the default and reproduces current selection.

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
