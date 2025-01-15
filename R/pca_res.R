#' Perform Principal Component Analysis (PCA) on Beta Values
#'
#' @description This function applies PCA to a matrix of beta values, visualizes the first four principal
#'  components (PC1 vs. PC2 and PC3 vs. PC4), and saves the resulting plots.
#'
#' @param beta_top100 A matrix of beta values. The rows represent probes, and the columns represent samples.
#' @param ss A data frame or data.table containing the sample sheet information.
#' @param sampGroups The name of the column in `ss` that contains sample group identifiers.
#' @param path A character string specifying the directory where plots will be saved (default: "./").
#' @param scale Logical, indicating whether to scale the variables (default: TRUE).
#' @param center Logical, indicating whether to center the variables (default: TRUE).
#' @param filename A prefix for the saved plot filenames (default: "")
#'
#' @return A `prcomp` object representing the results of PCA.
#'
#' @import ggplot2
#' @import ggfortify
#' @importFrom stats prcomp
#'
# @examples
# # Example usage:
# # pca_res(beta_top100, scale = TRUE, center = TRUE)
#'
#' @export
pca_res <- function(beta_top100,ss,sampGroups,path="./", scale = TRUE, center = TRUE,filename = "") {
  library(ggplot2)
  library(ggfortify)
   # Perform PCA on the transposed matrix
  pca_result <- prcomp(t(beta_top100), scale = scale, center = center)
  Sample_Group <- factor(ss[[sampGroups]])

  #Get proportion of variance explained bye each pca
  proportion<-summary(pca_result)$importance['Proportion of Variance',]

  #PC1 and PC2
  plt<-ggplot(pca_result,aes(x=PC1,y=PC2,color=ss[[sampGroups]]))+
    geom_point()+
    ggtitle('PC1 and PC2')+
    theme_bw()+ labs(color = 'Sample Group',
                          x = paste('PC1 (', proportion[1]*100, '%)', sep = ''),
                          y = paste('PC2 (', proportion[2]*100, '%)', sep = ''))
  save_plot(plt, path = path, filename = paste0(filename,'PC1PC2')) #"pcaplot_PC1PC2")

  #Plot PC3 and PC4
  plt34<-ggplot(pca_result,aes(x=PC3,y=PC4,color=ss[[sampGroups]]))+
    geom_point()+
    ggtitle('PC3 and PC4')+
    theme_bw()+ labs(color = 'Sample Group',
                     x = paste('PC3 (', proportion[3]*100, '%)', sep = ''),
                     y = paste('PC4 (', proportion[4]*100, '%)', sep = ''))
  save_plot(plt34, path = path,filename = paste0(filename,'PC3PC4'))


  return(pca_result)
}

