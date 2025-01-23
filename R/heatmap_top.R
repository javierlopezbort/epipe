#' Heatmap for beta values
#'
#' This function generates a heatmap of CpG probes using beta values.
#'
#' @param beta_values_matrix  A matrix of beta values (rows represent CpG probes, columns represent samples).
#' @param metadata A data frame containing sample metadata (e.g., sample group information).
#' @param sampGroups A string specifying the column name in `metadata` that represents the sample group for annotation.
#' @param pal Color palette for sample group (default: Brewer palette "Dark2")
#' @param path  A string indicating the directory path to save the generated heatmap plot. Default is the current directory. (default: "./")
#' @param filename Filename of the saved heatmap plot.
#'
#' @import ComplexHeatmap
#' @import RColorBrewer
#'
#' @return A heatmap object
#'
#' @export
#'
#' @examples
#' data("top_500_beta")
#' data("samplesheet")
#' heatmap_top(top_500_beta, samplesheet, sampGroups = "Condition")
#'
#'
heatmap_top<-function(beta_values_matrix,metadata,sampGroups,pal = RColorBrewer::brewer.pal(8, "Dark2"),
                      path = "./",filename = "heatmap"){


  # Add top annotations

  Condition <- factor(metadata[[sampGroups]])
  unique_conditions <- unique(Condition)
  color_palette <- pal[1:length(unique_conditions)]

  condition_colors<-setNames(color_palette,unique_conditions)

  ha <- HeatmapAnnotation(
    df = data.frame(Sample_Group = Condition),
    col = list(Sample_Group = condition_colors)
  )



  # Generate the heatmap
  plot<-Heatmap(beta_values_matrix,
          name='methylation',
          row_title = 'CpG probes',
          column_title='Samples',
          show_row_names = F,top_annotation = ha)

  # Save the plot
  save_plot(plot, path = path, filename = filename)

  return(plot)

}



