% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_category.R
\name{get_category}
\alias{get_category}
\title{Get Category from Input Sample Sheet}
\usage{
get_category(object)
}
\arguments{
\item{object}{A \code{data.table} or \code{data.frame} representing the sample sheet.}
}
\value{
The input object with an added \code{category} attribute that categorizes each column.
}
\description{
This function assigns a category to each variable of the sample sheet to ensure proper handling of
the information. Categories include: "ids" for identifiers, "batch" for batch variables, "covs" for covariates,
and "mgroups" for metadata groups.
}
\details{
\itemize{
\item \strong{Identifiers ("ids")}: Columns like \code{Sample_Name}, \code{barcode}, \code{Basename}, and \code{Sentrix_Position}.
\item \strong{Batch variables ("batch")}: Columns like \code{Sentrix_ID} and \code{batch}.
\item \strong{Covariates ("covs")}: Columns like \code{Type}, \code{Condition}, \code{Age}, \code{Sex}, \code{predictedSex}, and \code{predictedAge}.
\item \strong{Metadata groups ("mgroups")}: Columns like \code{Sample_Group}.
}

If a column cannot be categorized, it is labeled as "unknown" and a warning is issued. At least one column must belong to the "ids" category; otherwise, the function stops with an error.
If \code{Sentrix_ID} column is not found, it is extracted from \code{barcode}.
}
\examples{
samplesheet_path <- system.file("extdata", "EPICv2/samplesheet_EPICv2.rds", package = "epipe")
res <- readRDS(samplesheet_path)
sample_sheet <- get_category(res)
attributes(sample_sheet)

}
