# ADAPRE\_RECOMB demo

This folder contains a self-contained toy example of **ADAPRE** (adaptive row-wise penalties for inspre) together with the baseline inspre formulation, using a simple synthetic gene regulatory network.

## Repository layout

- `run_adapre_demo.R`: Main script to run the toy simulation and compare baseline inspre vs ADAPRE (VU/UV), including summary metrics and a degree-vs-true-degree plot.
- `R/adapre_inspre_core.R`: Core R implementation of inspre + helpers, extended with `adaptive_lambda` and instrument-strength-based penalties (ADAPRE).
- `R/simulate_two_interventions.R`: Data generator (from `generate_dataset_two_int_v2.R`) for two intervention-strength settings.
- `src/adapre_inspre_cpp.cpp`: C++ updates (via Rcpp) implementing the ADMM updates, including **both** `VU` and `UV` constraint paths with adaptive row-wise penalties.

## Requirements

R packages:

- `Rcpp`, `RcppEigen`
- `dplyr`, `purrr`, `ggplot2`
- `expm`, `igraph`, `nnls`, `caret` (needed by the underlying inspre helpers)
- `inspre` (for the network simulator used via `generate_dataset_two_int_v2.R`)

Install as needed, e.g.:

```r
install.packages(c("Rcpp", "RcppEigen", "dplyr", "purrr", "expm", "igraph", "nnls", "caret"))
devtools::install_github("brielin/inspre", INSTALL_opts = "--install-tests")
```

## Running the toy example

From inside this `ADAPRE_RECOMB` folder, run:

```bash
Rscript example_adapre_vs_inspre.R
```

The script will:

- Generate a small synthetic GRN with heterogeneous intervention strengths.
- Fit baseline inspre (`adaptive_lambda = FALSE`) and ADAPRE (`adaptive_lambda = TRUE`) under both `VU` and `UV` constraints.
- Print, for each setting, the best lambda (by training error) and the correlation between the true direct-effect matrix `G` and the estimated `R_hat` off-diagonal entries.
- Create and save a scatter plot (`degree_vs_true_degree_ADAPRE_VU.png`) showing, for each lambda, estimated vs true out-degree per gene, with points colored by intervention strength.
