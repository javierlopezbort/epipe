% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/apply_filter_dmrs.R
\name{make_ribbon_dt}
\alias{make_ribbon_dt}
\title{Make Ribbon Data Table}
\usage{
make_ribbon_dt(var, dt)
}
\arguments{
\item{var}{Variable for which ribbon data table is generated.}

\item{dt}{Data.table containing DMR information.}
}
\value{
A data.table for ribbon plotting.
}
\description{
This function generates a data.table for ribbon plotting based on a variable and data.table.
}
\examples{
# Example usage:
#ribbon_dt <- make_ribbon_dt(var = "min.cpg", dt = dt_data)

}
