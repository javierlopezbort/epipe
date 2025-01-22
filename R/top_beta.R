#' Extract Top Rows Based on Standard Deviation of Beta Values
#'
#' This function takes a matrix of beta values, computes the standard deviation
#' for each row, and extracts the top rows with the highest standard deviations.
#'
#' @param beta_values A matrix of beta values (rows represent CpG probes, columns represent samples).
#' @param bidx A list with two elements: `rn` (row names) and `cn` (column names) to set the row and column names for `beta_values`.
#' @param what A string specifying what to return. If "betas", returns the subset of rows with the highest standard deviations. If "rows", returns the row names (default is "betas").
#' @param n Number of top rows to extract (default: 1000).
#'
#' @return A matrix of beta values corresponding to the top `n` rows with the highest standard deviations
#'
#' @importFrom stats sd
#' @export
#'
#' @examples
#' # Example usage:
#' data("beta_matrix")
#' top_beta(beta_matrix, list(rn=rownames(beta_matrix),cn=colnames(beta_matrix)), n = 500)
#'



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

