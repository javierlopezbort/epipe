#' Density plot
#'
#' This function calculates beta values and generates a density plot.
#'
#' @param object An object containing intensity values (not betas).
#' @param metadata A data frame containing sample metadata.
#' @param sampGroups A character string specifying the column name in `metadata` that contains the sample grouping.
#' @param path A character string specifying the output directory for the density plot image. Defaults to "./".
#' @param norm_method A character string specifying the normalization method used.
#'
#' @return Generates a PNG file with the density plot in the specified path.
#'
#' @import minfi
#' @import grDevices
#' @export


denplot<-function(object, metadata, sampGroups, path="./",norm_method){
  object<-minfi::getBeta(object)
  grDevices::png(file = paste0(path, "density_plot_after_norm.png"))
  Sample_Group <- factor(metadata[[sampGroups]])
  title<-paste0('Beta values distribution after ',norm_method,' normalizationn')
  minfi::densityPlot(object, metadata[[sampGroups]],main = title)
  dev.off()
}

# Example
#denplot(filtered,ss,sampGroups,norm_method=norm_function)


