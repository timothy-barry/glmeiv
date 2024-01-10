
<!-- README.md is generated from README.Rmd. Please edit that file -->

# glmeiv

GLM-EIV is GLM-oriented variant of the errors-in-variables model.
GLM-EIV models the response as a GLM of an unobserved predictor. A noisy
version of the predictor is observed, which itself is modeled as a GLM
of the true predictor. This repository primarily is intended for reproduction of the results reported in the paper "Exponential
family measurement error models for single-cell CRISPR screens." Users who wish to employ the gRNA-only mixture assignment functionality
of GLM-EIV to assign gRNAs to cells should reference the [sceptre package](https://katsevich-lab.github.io/sceptre/).

# Installation

To install from Github use the following:

``` r
devtools::install_github("timothy-barry/glmeiv")
```
