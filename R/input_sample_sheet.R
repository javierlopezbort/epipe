#' @title Load data target factory.
#' @description Define 4 targets:
#' 1. Track the user-supplied data file.
#' 2. Read the data using `read_file()` .
#' 3. Generate category to handle info in sample_sheet.
#' @return A list of target objects.
#' @export
#' @param file Character, Sample sheet file path.
input_sample_sheet <- function(file) {
  list(
    targets::tar_target_raw("file", file, format = "file", deployment = "main"),
    targets::tar_target_raw("ss", quote(read_file(file)), deployment = "main"),
    targets::tar_target_raw("sample_sheet", quote(get_category(ss)),deployment = "main")   # Category attribute
  )
}
