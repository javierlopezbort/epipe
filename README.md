# EPIPE: (E)pigenetics (PIPE)line

**EPIPE** is a bioinformatics pipeline designed to analyze methylation data from microarrays. It automates the entire workflow, from processing raw data to identifying and annotating Differentially Methylated Positions (DMPs) and Regions (DMRs).

Key Features:

-   Fully automated pipeline for methylation analysis

-   Built on the targets framework for efficient and reproducible execution

-   Customizable configuration file for easy parameter adjustments

-   Users can modify and extend the pipeline for custom workflows

## Installation

You can install the development version of epipe with:

``` r
# install.packages("devtools")
options(timeout = 600)
devtools::install_github("ijcBIT/epipe",ref='development')
```

## Quick start

You can run the pipeline using the following commands:

``` r
library(epipe)
targets::tar_make()
```
