# isotwas 2.0.0.9000

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

## Improvements

* `compute_isotwas()` records failed methods in its result. `isotwas_result`
  gains a `failed_methods` data frame of the method name and the message it
  failed with, and `print()` reports it. Failures were previously visible only
  as a `verbose` console message.

* `isotwas_methods()` is the single source of truth for the method names, used
  both to build the comparison columns and to check each model's label.
