#' Rename RGSet and Adjust Column Names
#'
#' This function renames an RGSet object and adjusts column names based on target information.
#'
#' @param res An RGSet object to be processed.
#' @param targets Target information for column renaming.
#' @param newname New name for column renaming (optional).
#' @param exclude Columns to exclude from the result.
#' @param idcol Identifier column for matching with the RGSet columns.
#'
#' @return Modified RGSet object with adjusted column names and optional renaming.
#'
#' @import Biobase
#' @import SummarizedExperiment
#' @import data.table
#'
#' @examples
#' # Rename RGSet with targets information
#' res <- readRDS("inst/extdata/rgSet.rds")
#' targets <- readRDS("inst/sample_sheet.rds")
#' renamed_data <- name_rgset(res, targets, newname = "new_name", exclude = c("excluded_col1", "excluded_col2"))
#'
#' @export
name_rgset <- function(res, targets, newname = NULL, exclude = NULL, idcol = "barcode") {
  if (!requireNamespace("Biobase", quietly = TRUE) || !requireNamespace("SummarizedExperiment", quietly = TRUE) || !requireNamespace("data.table", quietly = TRUE)) {
    stop("Required packages 'Biobase', 'SummarizedExperiment', or 'data.table' are not installed.")
  }

  requireNamespace("Biobase")
  requireNamespace("SummarizedExperiment")
  requireNamespace("data.table")

  targets <- as.data.table(targets)
  if(is.null(idcol))idcol <- names(targets)[attributes(targets)$category=="ids"][1]
  cn <- targets[[idcol]]
  colnames(res) <- cn
  colnames(res@assays@data$Green) <- cn
  colnames(res@assays@data$Red) <- cn

  data.table::setkey(targets, "barcode")

  pheno <- as(targets, "DataFrame")
  rownames(pheno) <- cn

  stopifnot(rownames(pheno) == colnames(res))

  res@colData <- pheno

  if (!is.null(newname)){
    if(newname %in% colnames(pheno) & !any(duplicated(pheno[[newname]]))){
      colnames(res) <- res@colData[[newname]]
      rownames(pheno)<-colnames(res)
      res@colData <- pheno
    }

  }

  if (!is.null(exclude)) res <- res[, !colnames(res) %in% exclude]

  return(res)
}
