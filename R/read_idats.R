#' Read and Annotate IDAT Files
#'
#' Load IDAT files into a rgSet, names them correctly, and annotates them based on the array type.
#' @name read_idats
#' @param ss A data.table containing the sample sheet file.
#' @param arraytype A character string specifying the array platform used. (options: "450K", "EPIC", "EPICv2"). If `NULL`, the array type will be inferred from the sample sheet if available.
#' @param idats_folder Optional. A character string specifying the path to the folder containing IDAT files, if it is not correctly specified in the `Basename` column of the samplesheet.
#' @param idcol The column name in the sample sheet that contains unique identifiers for samples.(eg. Sample_Name)
#' @param description A description for the rgSet (default: "").
#' @param author The author of the rgSet (default: "ijcBIT").
#'
#' @return A SummarizedExperiment object.
#'
#' @details
#' The "Basename" column contains the path to the idat files (without _Red.idat and Grn.idat)
#' If idats_folder is provided, it will generate the "Basename" column by combining the `idats_folder`, `Sentrix_ID`,
#' and `Sentrix_Position` columns.
#' If idats_folder is not provided, the function reads the methylation data with the path that is in the `Basename` column.
#'
#'
#' @export
#' @import SummarizedExperiment
#' @importFrom data.table data.table
#' @import S4Vectors
#' @import dplyr
#'
#' @examples
#' \dontrun{
#'
#' # EXAMPLE 1: Basename column in the samplesheet
#' # Load the sample sheet
#' data("samplesheet")
#'
#' rgSet <- read_idats(ss = samplesheet,
#'                     arraytype = 'EPICv2',
#'                     idats_folder = NULL,
#'                     idcol = "Sample_Name")
#'
#'
#' # EXAMPLE 2: Idats_folder is provided
#'
#' data("samplesheet")
#'
#' # Define the path to the folder containing the IDAT files
#' idats_folder <-system.file("extdata/EPICv2/idats/", package = "epipe")
#'
#' rgSet <- read_idats(ss = samplesheet,
#'                     arraytype = 'EPICv2',
#'                     idats_folder = idats_folder,
#'                     idcol = "Sample_Name")
#' }


read_idats <- function(ss, arraytype = NULL, idats_folder = NULL, idcol = NULL, description = "", author = "ijcBIT") {
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


  if (!is.null(idats_folder)) {
    # If idats_folder is provided, construct Basename using it
    ss <- ss %>%
      mutate(Basename = paste(idats_folder, Sentrix_ID, '_', Sentrix_Position, sep = ''),
             barcode = paste(Sentrix_ID, '_', Sentrix_Position, sep = ''))
  } else if (!'Basename' %in% names(ss)) {
    # If idats_folder is not provided, check for Basename column
    stop("Neither 'idats_folder' nor 'Basename' column provided. Unable to construct Basename.")
  }

  rgSet <- read_meth(targets = ss, arraytype = arraytype)
  rgSet <- name_rgset(rgSet, ss, idcol = idcol) #Adds column names

  metadata(rgSet) <- list(
    description = description,
    author = author,
    category = unlist(attributes(ss)$category),
    date = Sys.Date()
  )
  return(rgSet)
}
