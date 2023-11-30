#' Create DataTable with custom options.
#'
#' This function creates a DataTable with custom options, including download buttons for the entire dataset.
#'
#' @param data A data frame to be displayed in the DataTable.
#'
#' @return A DataTable with custom options.
#'
#' @import DT
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' dtable(my_data)
#' }
custom_datatable <- function(data) {
  DT::datatable(
    { data },
    filter = 'top',
    fillContainer = FALSE,
    extensions = 'Buttons',
    options = list(
      paging = TRUE,
      pageLength = 10,
      searching = TRUE,
      fixedColumns = TRUE,
      autoWidth = FALSE,
      scrollX = TRUE,
      digits = 4,
      ordering = TRUE,
      dom = 'Bfrtip',
      buttons = list(
        list(
          extend = "collection",
          text = 'download entire dataset',
          action = DT::JS("function ( e, dt, node, config ) {
                           Shiny.setInputValue('test', true, {priority: 'event'});
                         }")
        ),
        'copy',
        'csv',
        'excel'
      ),
      class = "display",
      server = TRUE
    ),
  ) |> DT::formatRound(which(sapply(data, is.double)), 4)
}
