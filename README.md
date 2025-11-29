## Replication Files for [Synthetic Parallel Trends](https://drive.google.com/file/d/1MD1JSP1aNwMH1MtrSSLZH9HQFjY-bNlD/view?usp=sharing)
Please download the replication package, [synthdid-sdid-paper](https://www.openicpsr.org/openicpsr/project/146381/version/V1/view), provided by [Arkhangelsky, Athey, Hirshberg, Imbens, and Wager (2021)](https://www.aeaweb.org/articles?id=10.1257/aer.20190159). Place the entire downloaded folder `synthdid-sdid-paper` in the directory `code`.

In the directory `code`:
- `simulation-PT.R` replicates the results in Rows 1-2 of Table 1.
- `simulation-LFM.R` replicates the results in Rows 3-4 of Table 1.
- `all-func.R` provides the code for implementing the inference procedure in Section 4 and constructing the DGPs detailed in Section 5.
- `prep-data.R` produces `cps_rep_cross_sec.csv`, the cleaned CPS repeated cross-sectional sample used for generating data in the simulation designs.

All results will be directly printed to the RStudio console and also saved to the directory `result`. The results reported in Table 1 of the paper were run on a 2021 10-core M1 MacBook Pro with 16GB RAM under
```
R version 4.4.1 (2024-06-14)
Platform: aarch64-apple-darwin20
Running under: macOS Sonoma 14.6.1

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
```

## Note on Reproducibility
Even with an identical random R seed, exact numerical reproducibility across computing platforms is not guaranteed. Small discrepancies may arise from differences in floating-point arithmetic, BLAS/LAPACK implementations, package or R version differences, or the behavior of numerical solvers and multithreaded operations. 

For example, the first [arXiv](https://arxiv.org/abs/2511.05870v1) version of the paper used the default [OSQP](https://github.com/osqp) settings for optimization parameters that adaptively update based on computation time, which varies across platforms; see the discussion [here](https://github.com/osqp/osqp-r/issues/19#issuecomment-636954982) for more detail. The [latest version](https://drive.google.com/file/d/1MD1JSP1aNwMH1MtrSSLZH9HQFjY-bNlD/view?usp=sharing) of the paper fixes this source of variability by setting `adaptive_rho_interval=50` in `osqp::osqpSettings`, which disables the time-based adaptation. 

Another source of platform-specific variation appears in `simulation-PT.R` during the construction of DGP-2, which introduces a violation of post-treatment parallel trends. The relevant code snippet is:
```
for (trial in 1:20000){
  # try a different a_post that has more dispersion so we get a meaningful violation of post-trend
  a_post_alt <- a_post + rnorm(length(a_post), 0, 4*se_pt[length(se_pt)]) 
  # set violation of post-trend to be 2x se from event study
  get_w_alt <- get_alt_omega(A_pre, b_pre, a_post_alt, pi_0,
                             delta_target=2*se_pt[length(se_pt)])
  if(!is.null(get_w_alt$w_alt)) break
}
```
This loop searches for an alternative convex weight that (i) preserves parallel pre-trends and (ii) produces a post-treatment violation of parallel trends at least twice the DID standard error.
Whether `get_w_alt$w_alt` is `NULL` in a given `trial` depends on numerical factors such as BLAS/LAPACK precision (e.g., Apple Accelerate vs. OpenBLAS), small differences in QP solver tolerances, and finite-precision feasibility checks inside `get_alt_omega()`. If one platform exits the loop earlier than another (e.g., success at different values of `trial`), then that platform consumes a different number of `rnorm()` draws. This desynchronizes the global random number stream, so all subsequent random elements of that Monte Carlo replication follow different random sequences. 

The [latest version](https://drive.google.com/file/d/1MD1JSP1aNwMH1MtrSSLZH9HQFjY-bNlD/view?usp=sharing) of the paper reduces this source of platform-dependent randomness by setting a new seed `set.seed(iter)` at the beginning of each Monte Carlo iteration `iter`, preventing numerical differences inside the loop above from affecting the remainder of the simulation, and also increases the number of Monte Carlo replications from 500 to 1000 to further decrease the influence of platform-specific numerical differences. Any remaining discrepancies are expected to be minor and unlikely to affect the qualitative conclusions in Section 5.
