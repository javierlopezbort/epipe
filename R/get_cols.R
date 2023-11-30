#' Generate Categorical Color Palette
#'
#' Generates a categorical color palette based on the levels of a factor.
#'
#' @param factor_variable A factor variable.
#' @param pal Vector with colors (default is a set of categorical colors).
#' @return Vector with a color palette, same size as levels(factor_variable).
#' @keywords internal
#'
#' @examples
#' # Generate a color palette for a factor variable
#' factor_var <- factor(c("A", "B", "C", "A", "B"))
#' colors <- get_cols(factor_var)
#'
#' @import RColorBrewer
get_cols <- function(factor_variable, pal = NULL) {
  factor(factor_variable)
  levels_factor <- levels(factor_variable)

  if (is.null(pal)) {
    # Default color palette (modify as needed)
    if(length(levels_factor)<9){
      pal <- RColorBrewer::brewer.pal(length(levels_factor), "Dark2")
    }else{
      pal <- c(
        "#191919", "#0075DC", "#F0A0FF", "#993F00", "#005C31", "#FFE100", "#FF0010",
        "#FFCC99", "#808080", "#94FFB5", "#8F7C00", "#9DCC00", "#C20088", "#003380",
        "#FFA405", "#FFA8BB", "#426600", "#2BCE48", "#4C005C", "#00998F", "#E0FF66",
        "#740AFF", "#990000", "#FFFF80", "#5EF1F2", "#FF5005"
      )
    }
  }

  cols <- pal[as.numeric(factor_variable)]
  cols <- setNames(cols, labels(factor_variable))

  return(cols)
}
