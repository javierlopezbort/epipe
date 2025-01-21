#' Density plot
#'
#' This function generates a density plot.
#'
#' @param object An object containing normalized values (not betas, they are computed internally).
#' @param metadata A data frame containing sample metadata.
#' @param sampGroups A character string specifying the column name in `metadata` that contains the sample grouping.
#' @param path A character string specifying the output directory for the density plot image. Defaults to "./".
#' @param norm_method A character string specifying the normalization method used for normalization.
#'
#' @return Generates a PNG file with the density plot in the specified path.
#'
#' @import minfi
#' @import grDevices
#' @export
#'
#' @examples
#' data(samplesheet)
#' data(mSet_normalized)
#' denplot(mSet_normalized,samplesheet,sampGroups='Type',norm_method='swan')
#'


denplot<-function(object, metadata, sampGroups, path="./",norm_method){
  object<-minfi::getBeta(object)
  grDevices::png(file = paste0(path, "density_plot_after_norm.png"))
  Sample_Group <- factor(metadata[[sampGroups]])
  title<-paste0('Beta values distribution after ',norm_method,' normalizationn')
  minfi::densityPlot(object, metadata[[sampGroups]],main = title)
  dev.off()
}


