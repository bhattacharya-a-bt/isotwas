# Time each isoTWAS fitting method and compare it with a reference build.
#
# Usage, from the package root:
#
#   Rscript benchmarks/bench_methods.R ref <n> <p> <q> <methods> [seed] [outdir]
#   Rscript benchmarks/bench_methods.R new <n> <p> <q> <methods> [seed] [outdir]
#
# "ref" runs the isotwas installed in the library (the build to compare
# against). "new" runs this source tree through pkgload::load_all(). Run "ref"
# first: "new" reads the reference fits from <outdir> (default
# bench_out/) and reports, per method, the speedup and the largest
# differences in weights, cross-validated predictions and R-squared.
#
# If the environment variable ISOTWAS_LIB names a library, either mode loads
# isotwas from it instead. Use this to time two installed builds:
# load_all() compiles the C code without optimisation (-O0), which slows the
# compiled solvers.
#
# <methods> is a comma-separated subset of:
#   multi_enet joinet spls uni_enet uni_blup uni_susie sgl mtlasso
#   stacking graph_reg mrce iso_mrce iso_glmnet
# iso_mrce is compute_isotwas(method = 'mrce_lasso') with 10 Omegas, the way
# compute_isotwas runs MRCE; mrce is one compute_mrce() call. iso_glmnet is
# compute_isotwas() with the four glmnet-based methods (multi_enet, joinet,
# stacking, univariate), so it includes what the methods share within one run
# (folds and compressed designs); against a build without shared folds its
# weights and R-squared differ, since the folds differ.
#
# Data are simulated: binomial dosages with MAF ~ U(0.05, 0.5), 10 causal
# SNPs shared by all isoforms, and correlated isoform noise. The timings
# transfer; the R-squared values say nothing about real data.

args <- commandArgs(TRUE)
if (length(args) < 5) stop("usage: bench_methods.R <ref|new> <n> <p> <q> <methods> [seed] [outdir]")
mode <- args[1]
n <- as.integer(args[2]); p <- as.integer(args[3]); q <- as.integer(args[4])
meths <- strsplit(args[5], ",")[[1]]
dseed <- if (length(args) >= 6) as.integer(args[6]) else 1L
outdir <- if (length(args) >= 7) args[7] else "bench_out"
dir.create(file.path(outdir, mode), recursive = TRUE, showWarnings = FALSE)

if (!mode %in% c("ref", "new")) stop("mode must be 'ref' or 'new'")
lib <- Sys.getenv("ISOTWAS_LIB")
if (nzchar(lib) || mode == "ref") {
  lib_loc <- if (nzchar(lib)) lib else NULL
  suppressMessages(library(isotwas, lib.loc = lib_loc))
  desc <- utils::packageDescription("isotwas", lib.loc = lib_loc)
  message(mode, ": isotwas ", desc$Version,
          if (!is.null(desc$RemoteSha)) paste0(" (", desc$RemoteSha, ")"),
          " from ", dirname(dirname(system.file(package = "isotwas"))))
} else {
  suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE))
}
ns <- asNamespace("isotwas")
g <- function(f) get(f, envir = ns)
# ISOTWAS_GLMNET_COMPRESS=TRUE (or FALSE) sets options(isotwas.glmnet_compress);
# unset, the glmnet-based methods compress only inside compute_isotwas with
# shared folds, so the per-method cases need TRUE to time compression
compress <- Sys.getenv("ISOTWAS_GLMNET_COMPRESS")
if (nzchar(compress)) options(isotwas.glmnet_compress = as.logical(compress))

set.seed(dseed)
maf <- runif(p, .05, .5)
X <- sapply(maf, function(f) rbinom(n, 2, f)); colnames(X) <- paste0("rs", 1:p)
B <- matrix(0, p, q); B[sample(p, 10), ] <- rnorm(10 * q, 0, .3)
Y <- X %*% B + matrix(rnorm(n * q), n, q) %*% chol(0.5 * diag(q) + 0.5)
colnames(Y) <- paste0("tx", 1:q)
om <- g("compute_omega")(Y, Y, 1, seq_len(n), method = "mean", nlambda = 5, verbose = FALSE)
Om <- om$icov[[5]]
S <- matrix(0.3, q, q); diag(S) <- 1

f <- list(
  multi_enet = function() g("multivariate_elasticnet")(X, Y, Om, nfolds = 5, verbose = FALSE, seed = 42),
  joinet = function() g("multivariate_joinet")(X, Y, Om, nfolds = 5, verbose = FALSE, seed = 42),
  spls = function() g("multivariate_spls")(X, Y, nfolds = 5, verbose = FALSE, par = FALSE, n.cores = 1, seed = 42),
  uni_enet = function() g("univariate_elasticnet")(X, Y, Om, nfolds = 5, verbose = FALSE, seed = 42),
  uni_blup = function() g("univariate_blup")(X, Y, Om, nfolds = 5, verbose = FALSE, seed = 42),
  uni_susie = function() g("univariate_susie")(X, Y, Om, nfolds = 5, verbose = FALSE, seed = 42),
  sgl = function() g("multivariate_sgl")(X, Y, nfolds = 5, verbose = FALSE, seed = 42),
  mtlasso = function() g("multivariate_mtlasso")(X, Y, regularization = "L21", nfolds = 5, verbose = FALSE, seed = 42),
  stacking = function() g("multivariate_stacking")(X, Y, Y.rep = Y, R = 1, id = seq_len(n), Omega = Om,
                                                   nfolds_stack = 5, verbose = FALSE, seed = 42),
  graph_reg = function() g("multivariate_graph_reg")(X, Y, similarity_matrix = S, nfolds = 5, verbose = FALSE, seed = 42),
  iso_mrce = function() g("compute_isotwas")(X, Y, Y.rep = Y, R = 1, id = seq_len(n), omega_est = "mean",
                                             omega_nlambda = 10, method = "mrce_lasso", run_all = FALSE,
                                             verbose = FALSE, seed = 42)$best_models,
  iso_glmnet = function() g("compute_isotwas")(X, Y, Y.rep = Y, R = 1, id = seq_len(n), omega_est = "mean",
                                               omega_nlambda = 5,
                                               method = c("multi_enet", "joinet", "stacking", "univariate"),
                                               run_all = FALSE, verbose = FALSE, seed = 42)$best_models,
  mrce = function() g("compute_mrce")(X, Y, lambda = NULL, nlambda = 20, Omega = Om, nfolds = 5,
                                      tol.in = 1e-3, maxit.in = 1e3, verbose = FALSE, seed = 42,
                                      parallel = FALSE, n.cores = 1)
)
unknown <- setdiff(meths, names(f))
if (length(unknown)) stop("unknown method(s): ", paste(unknown, collapse = ", "))

compare <- function(a, b) {
  w <- function(m) lapply(m$transcripts, function(t) setNames(t$weights$Weight, t$weights$SNP))
  wd <- max(mapply(function(x, y) {
    s <- union(names(x), names(y))
    if (!length(s)) return(0)
    xx <- setNames(numeric(length(s)), s); yy <- xx
    xx[names(x)] <- x; yy[names(y)] <- y
    max(abs(xx - yy)) / max(1e-12, max(abs(yy)))
  }, w(a), w(b)))
  pd <- max(mapply(function(x, y)
    max(abs(x$predicted - y$predicted)) / max(1e-12, max(abs(y$predicted))),
    a$transcripts, b$transcripts))
  c(weights_relerr = wd, pred_relerr = pd,
    r2_absdiff = max(abs(a$summary$r2 - b$summary$r2)))
}

tag <- sprintf("%d_%d_%d_%d", n, p, q, dseed)
for (m in meths) {
  invisible(utils::capture.output(
    t <- system.time(r <- try(f[[m]](), silent = TRUE))[["elapsed"]]))
  saveRDS(list(time = t, res = r), file.path(outdir, mode, sprintf("%s_%s.rds", m, tag)))
  line <- sprintf("%-5s %-10s n=%d p=%d q=%d %9.2f s", mode, m, n, p, q, t)
  if (inherits(r, "try-error")) line <- paste(line, "ERROR:", conditionMessage(attr(r, "condition")))
  ref_file <- file.path(outdir, "ref", sprintf("%s_%s.rds", m, tag))
  if (mode == "new" && file.exists(ref_file) && !inherits(r, "try-error")) {
    ref <- readRDS(ref_file)
    if (!inherits(ref$res, "try-error")) {
      d <- compare(r, ref$res)
      line <- paste(line, sprintf("| %6.1fx | weights %.1e | pred %.1e | r2 %.1e",
                                  ref$time / t, d[1], d[2], d[3]))
    }
  }
  cat(line, "\n")
}
