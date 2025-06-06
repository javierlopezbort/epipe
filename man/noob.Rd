% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalization_methods.R
\name{noob}
\alias{noob}
\title{Apply Noob Normalization to an RGChannelSet Object}
\usage{
noob(rgSet)
}
\arguments{
\item{rgSet}{RGChannelSet object to normalize.}
}
\value{
Normalized RGChannelSet object.
}
\description{
This function applies the Noob normalization method to an RGChannelSet object.The Noob method (Negative Control
Probe Normalization) is used to preprocess methylation data, adjusting for dye bias in the data.
}
