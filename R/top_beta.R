#' Extract Top Rows Based on Standard Deviation of Beta Values
#'
#' This function takes a matrix of beta values, computes the standard deviation
#' for each row, and extracts the top rows with the highest standard deviations.
#'
#' @param beta_values A matrix of beta values.
#' @param bidx A data frame with row and column names for beta_values.
#' @param n Number of top rows to extract (default: 1000).
#' @return A subset of rows from beta_values with the highest standard deviations.
#'
#' @examples
#' # Example usage:
#' # top_beta(matrix_data, data.frame(rn = rownames(matrix_data), cn = colnames(matrix_data)), n = 500)
#'
#' @importFrom stats sd
#'
#' @export
top_beta <- function(beta_values, bidx,what="betas", n = 1000) {
  beta_values <- beta_values[]
  # Set row and column names
  rownames(beta_values) <- bidx$rn
  colnames(beta_values) <- bidx$cn

  # Compute standard deviation for each row
  sdv <- apply(beta_values, 1, sd)

  # Extract top n rows based on standard deviation
  top_n <- names(head(sort(sdv, decreasing = TRUE), n))
  beta_top_n <- beta_values[top_n, ]
  if(what == "betas"){
    return(beta_top_n)
  }
  else{
    return(top_n)
  }
}

