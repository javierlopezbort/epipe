#' Guess type of DNA methylation array.
#'
#' Will report the platform of a DNA methylation array given the number
#' of probes measured. Adapted from an internal function in \code{minfi} package.
#' Returns \emph{unknown} if array platform cannot be determined.
#'
#' @param n_probes Number of probes detected
#' @return Array type, one of:
#' \itemize{
#'   \item IlluminaHumanMethylationEPIC
#'   \item IlluminaHumanMethylation450k
#'   \item IlluminaHumanMethylation27k
#'   \item HorvathMammalMethylChip40
#'   \item unknown
#' }
#' @keywords internal

guess_arraytype <- function(n_probes) {
  if (n_probes >= 622000 && n_probes <= 623000) {
    array = "IlluminaHumanMethylation450k"
  } else if (n_probes >= 1050000 && n_probes <= 1053000) {
    # NOTE: "Current EPIC scan type"
    array = "IlluminaHumanMethylationEPIC"
  } else if (n_probes >= 1032000 && n_probes <= 1033000) {
    # NOTE: "Old EPIC scan type"
  } else if (n_probes >= 54000 && n_probes <= 56000) {
    array = "IlluminaHumanMethylation27k"
  } else if (n_probes >= 41000 & n_probes <= 41100) {
    array = "HorvathMammalMethylChip40"
  } else {
    array = "Unknown"
  }
  return(array)
}
