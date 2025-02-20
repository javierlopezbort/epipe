# EPIPE: (E)pigenetics (PIPE)line

**EPIPE** is a bioinformatics pipeline designed to analyze methylation data from microarrays. It automates the entire workflow, from processing raw data to identifying and annotating differentially methylated positions (DMPs) and regions (DMRs).

Based on the `targets` toolkit, EPIPE ensures efficient and reproducible execution of every analysis step. It includes a customizable configuration file for easy parameter adjustments, allowing users to modify and extend the pipeline to fit their specific analysis needs.

## Installation

Install the development version of epipe with:

``` r
# install.packages("devtools")
options(timeout = 600)
devtools::install_github("ijcBIT/epipe")
```

## Quick start

You can run the pipeline using the following commands:

``` r
library(epipe)
epipe::run()
```
