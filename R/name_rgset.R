#' Adjust Column Names
#'
#' This function adjusts column names based on target information.
#'
#' @param res An RGSet object to be processed.
#' @param targets Target information for column renaming.
#' @param idcol Identifier column for matching with the RGSet columns.
#'
#' @return Modified RGSet object with adjusted column names.
#'
#' @import Biobase
#' @import SummarizedExperiment
#' @import data.table
#'
# @examples
# # Rename RGSet with targets information
# res <- readRDS("inst/extdata/rgSet.rds")
# targets <- readRDS("inst/sample_sheet.rds")
# renamed_data <- name_rgset(res, targets, idcol='basename')
#

name_rgset <- function(res, targets, idcol = "Sample_Name") {
  if (!requireNamespace("Biobase", quietly = TRUE) || !requireNamespace("SummarizedExperiment", quietly = TRUE) || !requireNamespace("data.table", quietly = TRUE)) {
    stop("Required packages 'Biobase', 'SummarizedExperiment', or 'data.table' are not installed.")
  }

  targets <- as.data.table(targets)
  if(is.null(idcol))idcol <- names(targets)[attributes(targets)$category=="ids"][1]
  cn <- targets[[idcol]]
  colnames(res) <- cn
  colnames(res@assays@data$Green) <- cn
  colnames(res@assays@data$Red) <- cn

  #data.table::setkey(targets, "Sample_Name")

  pheno <- as(targets, "DataFrame")
  rownames(pheno) <- cn

  stopifnot(rownames(pheno) == colnames(res))

  res@colData <- pheno

  return(res)
}
