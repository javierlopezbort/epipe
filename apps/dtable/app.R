#' Shiny App with DataTable and Download Buttons.
#'
#' This Shiny app creates a DataTable with custom options, including download buttons for the entire dataset.
#'
#' @param data A data frame to be displayed in the DataTable.
#'
#' @import DT
#' @import shiny
#' @importFrom shiny downloadHandler observeEvent showModal
#' @importFrom writexl write_xlsx
#'
#' @export
#'
#' @return A Shiny app with a DataTable and download buttons.
#'
#' @examples
#' \dontrun{
#' dtable(my_data)
#' }
library(shiny)
library(DT)
library(writexl)


dtable <- function(data) {

    server <- function(input, output, session) {

        renderDT <- function(data) {
            output$dt <- DT::renderDataTable({
                custom_datatable(data)
            })
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
        }
    }

    ui <- fluidPage(
        title = 'Select Table Rows',
        h1('A Client-side Table'),
        fluidRow(
            column(12, DT::dataTableOutput('dt'))
            # column(6, plotOutput('x2', height = 500))
        )
    )
    shinyApp(ui, server)
}
