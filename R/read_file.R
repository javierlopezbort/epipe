#' Read File as RDS, Using fread, or readxl for Excel files
#'
#' This function attempts to open a file, first checking if it is an RDS file. If not,
#' it tries to open the file using `data.table::fread()` for CSV files or `readxl::read_excel()` for Excel files.
#'
#' @param file_path The path to the file to be read.
#' @return Data read from the file if successful.
#'
#' @importFrom data.table fread
#' @importFrom readxl read_excel
# @examples
# # Read a file as RDS, using fread, or readxl
# read_file("path/to/your/file")
#'
read_file <- function(file_path) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Required package 'data.table' is not installed.")
  }

  rds_data <- tryCatch({
    readRDS(file_path)
  }, error = function(e) {
    NULL
  })

  if (!is.null(rds_data)) {
    message("File opened as RDS.")
    return(rds_data)
  } else {
    fread_data <- tryCatch({
      data.table::fread(file_path)
    }, error = function(e) {
      NULL
    })

    if (!is.null(fread_data)) {
      message("File opened using fread.")
      return(fread_data)
    } else {
      excel_data <- tryCatch({
        readxl::read_excel(file_path)
      }, error = function(e) {
        NULL
      })

      if (!is.null(excel_data)) {
        message("File opened using readxl.")
        return(excel_data)
      } else {
        stop("Unable to open the file as RDS, using fread, or readxl. Please use a valid file format.")
      }
    }
  }
}

