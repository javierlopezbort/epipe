#' Render DataTable and handle download events for internal use.
#'
#' This function renders a DataTable and handles download events for CSV and Excel formats.
#' It is intended for internal use within the package and should not be exported to the end user.
#' @name dtable
#' @title Dynamic table with download buttons
#' @param data A data frame to be rendered in the DataTable.
#'
#' @keywords internal
#'
#' @import DT
#' @import shiny
#' @importFrom shiny renderDataTable observeEvent showModal downloadHandler
#' @importFrom writexl write_xlsx
#'
#' @return A DataTable with formatted data and handlers for CSV and Excel downloads.
#'
#' @examples
#' \dontrun{
#' dtable(my_data)
#' }
library(DT)
library(writexl)
library(shiny)

dtable <- function(data) {
  mt <- DT::renderDataTable(
    custom_datatable(data) |> formatRound(which(sapply(data, is.double)), 4)
  )

  shiny::observeEvent(input$test, {
    print("hello")
    showModal(myModal())
  })

  output$download1 <- shiny::downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(data, file)
    }
  )

  output$download2 <- shiny::downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".xlsx", sep="")
    },
    content = function(file) {
      writexl::write_xlsx(data, file)
    }
  )

  mt
}
