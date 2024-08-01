#' Boxplot and statistical test
#' 
#' Performs t.test to check for statistical differences for cell type deconvolution between two groups and creates a boxplot.
#' 
#' @title ttest_boxplot
#' @param data Data frame containing metadata
#' @param sampGroup Name of the variable in the metadata that refers to the sample Grops or condition.
#' @param path Path to save the plots (default: "./")
#' 
#' 
#' @export

ttest_boxplot<-function(data,sampGroup,path = "./",filename='Cell_type_prop.png'){
  
  #libraries needed: 
  library(tidyr)
  library(ggpubr)
  library(rstatix)
  
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
  bxp <- ggboxplot(
    df_long, x = "cell_type", y = "value",
    fill = "condition", palette = c("#00AFBB", "#FF7777")
  )
  
  # Add statistical test results to the plot
  plot <- bxp + stat_pvalue_manual(
    stat.test, label = "{p} {p.signif}", tip.length = 0.01, label.vjust = -0.1
  ) + labs(title = 'Boxplot. Cell type proportions')
 
  # Save the plot
  file_path=paste0(path,filename)
  ggsave(file_path, plot, width = 10, height = 7)
   
}

# Use: 

#ttest_boxplot(ss_clean_allvariables_ex_quantile_sex,sampGroup = 'condition')
