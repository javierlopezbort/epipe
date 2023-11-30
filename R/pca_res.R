#' Perform Principal Component Analysis (PCA) on Beta Values
#'
#' This function applies Principal Component Analysis (PCA) to a matrix of beta values.
#'
#' @param beta_top100 A matrix of beta values.
#' @param scale Logical, indicating whether to scale the variables (default: TRUE).
#' @param center Logical, indicating whether to center the variables (default: TRUE).
#' @return A prcomp object representing the results of PCA.
#'
#' @importFrom stats prcomp
#'
#' @examples
#' # Example usage:
#' # pca_res(beta_top100, scale = TRUE, center = TRUE)
#'
#' @export
pca_res <- function(beta_top100, scale = TRUE, center = TRUE) {
  # Perform PCA on the transposed matrix
  pca_result <- prcomp(t(beta_top100), scale = scale, center = center)

  return(pca_result)
}
