% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/apply_filter_dmps.R
\name{filter_dmps}
\alias{filter_dmps}
\title{Filter Differentially Methylated Positions (DMPs)}
\usage{
filter_dmps(dmps, adj.p.value = 0.05, p.value = 0.01, mDiff = 0.3, s = F)
}
\arguments{
\item{dmps}{Data table containing DMPs.}

\item{adj.p.value}{Threshold for the adjusted p-value (default: 0.05).}

\item{p.value}{Threshold for the raw p-value (default: 0.01).}

\item{mDiff}{Threshold for the absolute mean methylation difference (default: 0.3).}

\item{s}{Logical indicating whether to include summary statistics.}
}
\value{
A data table of filtered DMPs
}
\description{
This function filters DMPs based on adjusted p-value, p-value, and mean methylation difference thresholds.
}
