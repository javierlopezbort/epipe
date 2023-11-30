#' Read IDAT Files and Assign Array Type
#'
#' Reads IDAT files and assigns the array type based on the provided information or guessed from the data.
#'
#' @inheritParams read_idats
#' @return RGSet object with assigned array type.
#' @export
#' @import minfi
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
