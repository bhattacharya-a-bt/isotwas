# Method benchmarks

`bench_methods.R` times each fitting method on one simulated gene, once with
a reference build of isotwas and once with a new one, and reports the
speedup and the largest differences between the two fits (weights,
cross-validated predictions, R-squared). This directory is not part of the
installed package.

From the package root, with the reference build installed:

```
Rscript benchmarks/bench_methods.R ref <n> <p> <q> <methods> [seed] [outdir]
Rscript benchmarks/bench_methods.R new <n> <p> <q> <methods> [seed] [outdir]
```

Run `ref` first; `new` reads its fits from `<outdir>` (default `bench_out/`)
and prints one comparison line per method. `<methods>` is a comma-separated
subset of `multi_enet joinet spls uni_enet uni_blup uni_susie sgl mtlasso
stacking graph_reg mrce iso_mrce iso_glmnet`.

By default `ref` loads the installed isotwas and `new` loads this source tree
with `pkgload::load_all()`. `load_all()` compiles the C code at `-O0`, which
slows the compiled solvers, so for timings install both builds and point each
run at its library with the environment variable `ISOTWAS_LIB`:

```
ISOTWAS_LIB=lib_ref Rscript benchmarks/bench_methods.R ref 1000 500 5 sgl,mrce
ISOTWAS_LIB=lib_new Rscript benchmarks/bench_methods.R new 1000 500 5 sgl,mrce
```

The glmnet-based methods compress their inputs by default only inside
`compute_isotwas()` with shared folds (`iso_glmnet`). To time compression for
the methods called on their own, set `ISOTWAS_GLMNET_COMPRESS=TRUE` for the
`new` run (it sets `options(isotwas.glmnet_compress)`).

Results for the speed-ups this script was written for are in
https://github.com/bhattacharya-a-bt/isotwas/pull/11.
