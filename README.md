# EPIPE: (E)pigenetics (PIPE)line

**EPIPE** is a bioinformatics pipeline designed to analyze methylation data from microarrays. It automates the entire workflow, from processing raw data to identifying and annotating Differentially Methylated Positions (DMPs) and Regions (DMRs).

Based on the `targets` framework, EPIPE ensures efficient and reproducible execution of every analysis step.

## Installation

You can install the development version of epipe with:

``` r
# install.packages("devtools")
devtools::install_github("ijcBIT/epipe",ref='development')
```

## Quick start

You can run the pipeline using the following commands:

``` r
library(epipe)
targets::tar_make()
```
