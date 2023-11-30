
#' Read and Annotate IDAT Files
#'
#' Loads IDAT files into an rgSet, names them correctly, and annotates them based on the array type.
#'
#' @param ss A data.table containing the sample sheet file.
#' @param arraytype Which platform has been used for the array c("450K","EPIC","EPICv2")
#' @return A list of target objects.
#' @export
#' @importFrom targets tar_target_raw
#' @importFrom methods as
read_idats <- function(ss, arraytype = NULL) {
  if (is.null(arraytype)) {
    if ("arraytype" %in% colnames(ss)) {
      if (length(unique(ss$arraytype)) > 1) {
        stop("arraytype must be consistent in all samples of rgSet object. Try subsetting rgSet by arraytype")
      } else arraytype <- unique(ss$arraytype)
    } else {
      arraytype <- "guess"
    }
  }

  arraytype <- switch(arraytype,
                      "EPICv2" = "IlluminaHumanMethylationEPICv2",
                      "EPIC" = "IlluminaHumanMethylationEPIC",
                      "450K" = "IlluminaHumanMethylation450k",
                      "27K" = "IlluminaHumanMethylation27k",
                      "mammal" = "HorvathMammalMethylChip40",
                      "guess")

  sym_arraytype <- as.symbol(arraytype)
  command_idats <- substitute(read_meth(targets = ss, arraytype = arraytype), env = list(arraytype = sym_arraytype))

  list(
    targets::tar_target_raw("nrgSet", command_idats),
    targets::tar_target_raw("rgSet", quote(name_rgset(nrgSet, ss, exclude = NULL, newname = idcol)),
                            deployment = "worker",
                            memory = "persistent")
  )
}


#' Read IDAT Files and Assign Array Type
#'
#' Reads IDAT files and assigns the array type based on the provided information or guessed from the data.
#'
#' @inheritParams read_idats
#' @return RGSet object with assigned array type.
#' @export
#' @import minfi
read_meth <- function(ss, arraytype = NULL) {
  # Read IDAT files
  rgSet <- minfi::read.metharray.exp(base = NULL, targets = ss, extended = TRUE, recursive = FALSE, verbose = FALSE, force = TRUE)

  # Guess array type if not provided or set to 'guess'
  if (is.null(arraytype) || arraytype == "guess") {
    arraytype <- guess_arraytype(nrow(rgSet))  # Replace guess_arraytype with your guessing function
  }

  # Assign array-specific annotation
  if (arraytype == "IlluminaHumanMethylation450k") {
    rgSet@annotation <- c(array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19")
  } else if (arraytype == "IlluminaHumanMethylationEPIC") {
    rgSet@annotation <- c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b2.hg19")
  } else if (arraytype == "IlluminaHumanMethylationEPICv2") {
    rgSet@annotation <- c(array = "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38")
  } else if (arraytype == "IlluminaHumanMethylation27k") {
    rgSet@annotation <- c(array = "IlluminaHumanMethylation27k", annotation = "ilmn12.hg19")
  } else if (arraytype == "HorvathMammalMethylChip40") {
    rgSet@annotation <- c(array = "HorvathMammalMethylChip40", annotation = "test.unknown")
  } else if (arraytype == "unknown") {
    rgSet@annotation <- c(array = "chip.unknown", annotation = "test.unknown")
  }

  return(rgSet)
}

