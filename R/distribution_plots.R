#' Distribution plots
#'
#' This function generates either boxplots or barplots for all variables in a data frame,
#' colored by the variable of interest (such as sample groups or conditions).
#'
#' @param data Data frame containing metadata
#' @param variable_interest Name of the variable in the metadata that refers to the sample groups or condition.
#' @param path Path to save the plots (default: "./")
#'
#' @import ggplot2
#' @export
#' @examples
#' \dontrun{
#' data <- data.frame(
#'   Group = factor(c("A", "B", "A", "B", "A", "B")),
#'   Age = c(30, 25, 40, 35, 28, 45),
#'   Gender = factor(c("M", "F", "M", "F", "M", "F"))
#' )
#' distribution_plots(data, variable_interest = "Group")
#'}

distribution_plots<-function(data,variable_interest,path = "./"){

  if ("Basename" %in% colnames(data)) {
    data <- subset(data, select = -Basename)
  }
  if ("barcode" %in% colnames(data)) {
    data <- subset(data, select = -barcode)
  }
  if ("yMed" %in% colnames(data)) {
    data <- subset(data, select = -yMed)
  }
  if ("xMed" %in% colnames(data)) {
    data <- subset(data, select = -xMed)
  }
  if ("Sample_Name" %in% colnames(data)) {
    data <- subset(data, select = -Sample_Name)
  }
  if ("Sentrix_position" %in% colnames(data)) {
    data <- subset(data, select = -Sentrix_position)
  }


  for (variable in names(data)){
    if (variable == variable_interest) next
    if (is.numeric(data[[variable]])) {
      # Create a boxplot for numerical variables
      p <- ggplot(data, aes_string(x = variable_interest, y = variable, fill = variable_interest)) +
        geom_boxplot() +
        labs(title = paste("Boxplot of", variable, "colored by", variable_interest))

      filename<-paste0('Boxplot_',variable,'_colored_by_',variable_interest,'.png')
    } else if (is.factor(data[[variable]]) || is.character(data[[variable]])) {
      # Create a barplot for categorical variables
      p <- ggplot(data, aes_string(x = variable, fill = variable_interest)) +
        geom_bar() +
        labs(title = paste("Barplot of", variable, "colored by", variable_interest))

      filename<-paste0('Barplot_',variable,'_colored_by_',variable_interest,'.png')
    } else {
      next
    }

    # Save the plot
    file_path=paste0(path,filename)
    ggsave(file_path, p, width = 10, height = 7)
  }

}
