#' Use Betas Object in Parallel and Save to Disk
#'
#' This function takes a beta matrix from a rgSet object and converts it into a
#' bigstatsr::FBM (Filebacked Big Matrix) for parallel processing. This can help
#' in efficient storage and manipulation of beta values, especially when the
#' matrix is large. The betas are saved to disk using the bigstatsr package.
#'
#' @param rgSet A SummarizedExperiment object containing beta values.
#' @return A bigstatsr::FBM object representing the beta matrix.
#'
#' @importFrom bigstatsr FBM
#' @importFrom minfi getBeta
#'
#' @export
#'
betas_disk <- function(rgSet, ...) {

  # Extract beta matrix from rgSet
  beta_normalized <- minfi::getBeta(rgSet)

  # Create a bigstatsr::FBM and save betas to disk
  betas <- bigstatsr::FBM(NROW(beta_normalized),NCOL(beta_normalized), ...);betas[] <- as.matrix(beta_normalized)
  # betas <- bigstatsr::FBM(nrow(beta_normalized), ncol(beta_normalized))
  # betas[] <- as.matrix(beta_normalized)

  return(betas)
}
