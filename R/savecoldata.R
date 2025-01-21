#' Save Column Data from a SummarizedExperiment Object
#'
#' This function extracts the column data from a SummarizedExperiment object,
#' drops unused levels from factors in the resulting data frame, and saves the
#' data frame to a CSV file. Optionally, it performs an assertion check to ensure
#' consistency between row names of the SummarizedExperiment object and column
#' names of the data frame.
#'
#' @param rgset A SummarizedExperiment object.
#' @param file The name of the file to save the data (default: "results/ss_clean.csv").
#' @param quote A logical value indicating whether to quote the output file (default: FALSE).
#' @param check_assertion Logical, perform an assertion check (default: TRUE).
#' @param sep The field separator for the output file (default: ",").
#' @param dir The directory where the file should be saved (default: ".").
#' @param verbose Logical, print informative messages (default: FALSE).
#'
#' @return A data frame containing the cleaned column data.
#'
#' @examples
#' \dontrun{
#' data("mSet_normalized")
#' savecoldata(mSet_normalized)
#' }
#'
#' @import S4Vectors
#' @import assertthat
#' @importFrom utils write.table
#'
#' @export
savecoldata <- function(rgset, file = "results/ss_clean.csv", quote = FALSE,
                        check_assertion = TRUE, sep = ",", dir = ".", verbose = FALSE) {
  # Check input class
  if (!inherits(rgset, "SummarizedExperiment")) {
    stop("Input must be a SummarizedExperiment object.")
  }

  # Extract and drop levels
  ss <- droplevels.data.frame(colData(rgset))

  # Assertion check
  if (check_assertion) {
    if(!assertthat::are_equal(colnames(rgset), rownames(ss)))stop("Assertion check failed.")
  }

  # Save file
  output_path <- file.path(dir, file)
  utils::write.table(ss, file = output_path, sep = sep, quote = quote, row.names = FALSE)

  # Informative message
  if (verbose) {
    cat("Column data saved successfully to:", output_path, "\n")
  }

  # Return the data frame
  return(ss)
}
