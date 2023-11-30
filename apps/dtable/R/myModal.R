#' Modal dialog with download buttons for internal use.
#'
#' This function creates a modal dialog with download buttons for CSV and Excel files.
#' It is intended for internal use within the package and should not be exported to the end user.
#'
#' @keywords internal
#'
#' @import shiny
#' @importFrom shiny modalDialog
#' @importFrom shiny downloadButton
#'
#' @return A modal dialog with download buttons.
#'
#' @export
#' @examples
#' \dontrun{
#' myModal()
#' }
myModal <- function() {

  shiny::modalDialog(
    shiny::downloadButton("download1", "Download data as csv"),
    br(),
    br(),
    shiny::downloadButton("download2", "Download data as excel"),
    easyClose = TRUE, title = "Download Table"
  )
}
