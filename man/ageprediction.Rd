% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ageprediction.R
\name{ageprediction}
\alias{ageprediction}
\title{Age prediction}
\usage{
ageprediction(object, clock_name = "Horvath", predict_age = TRUE)
}
\arguments{
\item{object}{MethylSet object}

\item{clock_name}{Name of the age clock used to predict age. Default=Horvath (horvath 2013). Currently "Hannum", "Levine", "BNN", "skinHorvath", "PedBE", "Wu", "TL", "BLUP", "EN" and "all" are available.}

\item{predict_age}{Logical value, whether to predict age or not. Default=TRUE}
}
\value{
mSet object with predicted age values in colData().
}
\description{
This function predicts biological age from DNA methylation data.
}
\examples{
data("mSet_normalized")
ageprediction(mSet_normalized,clock_name='Horvath')
}
