#' Read files
#'
#' @param file_path A string specifying the path to the file  (include it. eg. 'data/file.csv').
#' @return A data frame or R object, depending on the file type.
#' @importFrom tools file_ext
#'
#' @export

# Example usage:
# df <- read_file("data.csv")
# df <- read_file("data.rds")


read_file <- function(file_path) {

  #Determine the file extension
  file_ext <- tools::file_ext(file_path)

  # Read the file based on its extension
  if (file_ext == "rds") {
    return(readRDS(file_path))
  } else if (file_ext == "csv") {
    return(read.csv(file_path, stringsAsFactors = FALSE))
  } else {
    stop("Unsupported file type: ", file_ext)
  }

}
