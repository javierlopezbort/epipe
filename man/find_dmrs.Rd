% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_dmrs2.R
\name{find_dmrs}
\alias{find_dmrs}
\title{Find DMRs}
\usage{
find_dmrs(
  beta_values,
  model,
  fdr = 0.01,
  p.value = "fdr",
  bcutoff = 0.05,
  min.cpg = 5,
  pal,
  output_dir = "./",
  arraytype
)
}
\arguments{
\item{beta_values}{Beta values matrix.}

\item{model}{Linear model used for differential methylation analysis.}

\item{fdr}{False Discovery Rate threshold.}

\item{p.value}{P-value threshold (not recommended, use fdr).}

\item{bcutoff}{Beta value cutoff.}

\item{min.cpg}{Minimum number of CpGs in a DMR.}

\item{output_dir}{Path to folder to save results.}

\item{arraytype}{Type of array. EPICv2, EPIC or 450K.}
}
\value{
A data.table containing the identified DMRs.
}
\description{
Find Differentially Methylated Regions (DMRs) -2
}
\details{
This function identifies differentially methylated regions (DMRs) from beta values and a linear model.
}
