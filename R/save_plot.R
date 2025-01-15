#' Save Plot to File
#'
#' This function saves a plot to a file in PNG format.
#'
#' @param object A plot object (base R plot or ggplot2 plot).
#' @param filename Name of the file to save the plot.
#' @param width Width of the plot in inches (default: 480).
#' @param height Height of the plot in inches (default: 620).
#' @param path Directory path for saving the plot (default: "./").
#' @return The plot object.
#'
save_plot <- function(object, filename, width = 480, height = 620, path = "./") {
  # Create the directory if it doesn't exist
  dir.create(path, recursive = TRUE, showWarnings = FALSE)

  # Determine the type of plot (base R or ggplot2)
  is_ggplot <- inherits(object, "gg") || inherits(object, "ggplot")

  # Save the plot based on its type
  if (is_ggplot) {
    # For ggplot2 plots, use ggsave function
    ggsave(
      file.path(path, paste0(filename, ".png")),
      plot = object,
      width = width / 100,
      height = height / 100
    )
  } else {
    # For base R plots, use png and print
    grDevices::png(
      file = file.path(path, paste0(filename, ".png")),
      width = width,
      height = height,
      units = "px",
      res = 100
    )
    print(object)
    dev.off()
  }

  return(object)
}

