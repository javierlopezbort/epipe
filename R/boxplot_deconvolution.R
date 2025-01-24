#' Boxplot and statistical test
#'
#' This function performs a t-test to check for statistical differences in cell type proportions between two groups
#' and creates a boxplot with statistical annotations.
#'
#' @param data A data frame containing metadata and cell type proportions.
#' @param sampGroup A character string specifying the name of the variable in the metadata that refers to the sample groups or condition.
#' @param pal Color palette to use
#' @param path A character string specifying the directory to save the plot (default: "./").
#' @param filename A character string specifying the file name for the saved plot (default: "Cell_type_prop.png").
#'
#' @return Saves a boxplot as a PNG file.
#'
#' @import tidyr
#' @import ggpubr
#' @import rstatix
#'
#' @export
#'
#' @examples
#' data("ss_all_variables")
#' ttest_boxplot(ss_all_variables,sampGroup = 'Sample_Group')

ttest_boxplot<-function(data,sampGroup,pal=NULL,path = "./",filename='Cell_type_prop.png'){

  column_names <- c(sampGroup,"CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu", "leuk.prop")
  selected_data <- data[, column_names]

  colnames(selected_data)[colnames(selected_data) == sampGroup] <- "condition"

  #Convert it to a data frame
  df <- as.data.frame(selected_data)

  # Pivot the dataframe to long format
  df_long <- pivot_longer(df,
                          cols = -condition,  # Assuming 'condition' is the column name for comparison
                          names_to = "cell_type",
                          values_to = "value")

  # Convert necessary columns to factors
  df_long<-as.data.frame(df_long)
  df_long$cell_type<-as.factor(df_long$cell_type)
  df_long$condition<-as.factor(df_long$condition)

  # Perform statistical tests
  stat.test <- df_long %>%
    group_by(cell_type) %>%
    t_test(value ~ condition) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p")

  stat.test <- stat.test %>%
    add_xy_position(x = 'cell_type', dodge = 0.8)

  # Create the boxplot
  if (is.null(pal)) {
    # Default palette with 9 distinct colors
    pal <- c("#1B9E77", "#7570B3", "#D95F02", "#E7298A", "#66A61E",
                 "#E6AB02", "#A6761D", "#666666", "#FF0010")
  }

  palette_colors <- pal[seq_along(unique(df_long$condition))]

  bxp <- ggboxplot(
    df_long, x = "cell_type", y = "value",
    fill = "condition", palette = palette_colors
  )

  # Add statistical test results to the plot
  plot <- bxp + stat_pvalue_manual(
    stat.test, label = "{p} {p.signif}", tip.length = 0.01, label.vjust = -0.1
  ) + labs(title = 'Boxplot. Cell type proportions')

  # Save the plot
  file_path=paste0(path,filename)
  ggsave(file_path, plot, width = 10, height = 7)

}
