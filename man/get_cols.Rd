% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_cols.R
\name{get_cols}
\alias{get_cols}
\title{Generate Categorical Color Palette}
\usage{
get_cols(factor_variable, pal = NULL)
}
\arguments{
\item{factor_variable}{A factor variable.}

\item{pal}{Vector with colors (default is a set of categorical colors).}
}
\value{
Vector with a color palette, same size as levels(factor_variable).
}
\description{
Generates a categorical color palette based on the levels of a factor.
}
\keyword{internal}
