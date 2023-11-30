#' Read and Annotate IDAT Files
#'
#' Loads IDAT files into an rgSet, names them correctly, and annotates them based on the array type.
#'
#' @param ss A data.table containing the sample sheet file.
#' @param arraytype Which platform has been used for the array (options: "450K", "EPIC", "EPICv2").
#' @param idcol The column name in the sample sheet that contains unique identifiers for samples.
#' @param description A description for the rgSet (default: "").
#' @param author The author of the rgSet (default: "ijcBIT").
#'
#' @return A SummarizedExperiment object.
#'
#' @export
#' @import SummarizedExperiment
#' @importFrom data.table data.table

requireNamespace("S4Vectors")

read_idats <- function(ss, arraytype = NULL, idcol = NULL, description = "", author = "ijcBIT") {
  if (is.null(arraytype)) {
    if ("arraytype" %in% colnames(ss)) {
      if (length(unique(ss$arraytype)) > 1) {
        stop("arraytype must be consistent in all samples of rgSet object. Try subsetting rgSet by arraytype")
      } else arraytype <- unique(ss$arraytype)
    } else {
      arraytype <- "guess"
    }
  }

  arraytype <- switch(
    arraytype,
    "EPICv2" = "IlluminaHumanMethylationEPICv2",
    "EPIC" = "IlluminaHumanMethylationEPIC",
    "450K" = "IlluminaHumanMethylation450k",
    "27K" = "IlluminaHumanMethylation27k",
    "mammal" = "HorvathMammalMethylChip40",
    "guess"
  )
  rgSet <- read_meth(targets = ss, arraytype = arraytype)
  rgSet <- name_rgset(rgSet, ss, exclude = NULL, idcol = idcol)

  metadata(rgSet) <- list(
    description = description,
    author = author,
    category = unlist(attributes(ss)$category),
    date = Sys.Date()
  )
  return(rgSet)
}
