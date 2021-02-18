<!-- badges: start -->
[![R build status](https://github.com/HectorRDB/Ecume/workflows/R-CMD-check/badge.svg)](https://github.com/HectorRDB/Ecume/actions)
[![Codecov test coverage](https://codecov.io/gh/HectorRDB/Ecume/branch/main/graph/badge.svg)](https://codecov.io/gh/HectorRDB/Ecume?branch=main)
<!-- badges: end -->
  
# R package: Ecume

  
This package implemements (or re-implements in `R`) a variety of statistical tools used to do non-parametric two-sample (or k-sample) distribution comparisons in the univariate or multivariate case.

If you want to suggest a new test or a variant, open an [issue](https://github.com/HectorRDB/Ecume/issues).

## Installation

To install the current version of *Ecume*, run.

```
if(!requireNamespace("BiocManager", quietly = TRUE)) {
 install.packages("BiocManager") 
}
BiocManager::install("Ecume")
```

To install the development version in `R`, run 

```r
devtools::install_github("HectoRDB/Ecume")
```
