#' Read IDAT Files and Assign Array Type
#'
#' Reads IDAT files and assigns the array type based on the provided information or guessed from the data.
#'
#' @param targets A data frame or data.table containing sample metadata, including the paths to the IDAT files.
#' @param arraytype A character string specifying the array platform. Options include "450K", "EPIC", "EPICv2", "27K", or "mammal". If `NULL` or "guess", the array type will be inferred from the data.
#'
#' @return RGSet object with assigned array type.
#' @import minfi
#'
read_meth <- function(targets, arraytype = NULL) {
  # Read IDAT files
  rgSet <- minfi::read.metharray.exp(base = NULL, targets = targets, extended = TRUE, recursive = FALSE, verbose = FALSE, force = TRUE)

  # Guess array type if not provided or set to 'guess'
  if (is.null(arraytype) || arraytype == "guess") {
    arraytype <- guess_arraytype(nrow(rgSet))  # Replace guess_arraytype with your guessing function
  }
  rgSet <- set_array_annotation(rgSet, arraytype = arraytype)

  return(rgSet)
}
