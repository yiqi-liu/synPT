# Replication files for [Synthetic Parallel Trends](https://drive.google.com/file/d/1MD1JSP1aNwMH1MtrSSLZH9HQFjY-bNlD/view?usp=sharing)
Please download the replication package, [synthdid-sdid-paper](https://www.openicpsr.org/openicpsr/project/146381/version/V1/view), provided by [Arkhangelsky, Athey, Hirshberg, Imbens, and Wager (2021)](https://www.aeaweb.org/articles?id=10.1257/aer.20190159). Place the entire downloaded folder `synthdid-sdid-paper` in the directory `code`.

In the directory `code`:
- `simulation-PT.R` replicates the results in Rows 1-2 of Table 1.
- `simulation-LFM.R` replicates the results in Rows 3-4 of Table 1.
- `all-func.R` provides the code for implementing the inference procedure in Section 4 and constructing the DGPs detailed in Section 5.
- `prep-data.R` produces `cps_rep_cross_sec.csv`, the cleaned CPS repeated cross-sectional sample used for generating data in the simulation designs.

All results will be directly printed in the RStudio console and also saved to the directory `result`. The results in Table 1 of the paper is run on
```
R version 4.4.1 (2024-06-14)
Platform: aarch64-apple-darwin20
Running under: macOS Sonoma 14.6.1

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
```
