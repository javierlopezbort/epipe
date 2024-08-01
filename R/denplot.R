#' Density plot after normalization
#' @return dplot
#' @import Minfi
#' @export


denplot<-function(object, metadata, sampGroups, path="./",norm_method){
  library(minfi)
  object<-minfi::getBeta(object) # sharuia de posar: if class(object)==.. {}
  grDevices::png(file = paste0(path, "density_plot_after_norm.png"))
  Sample_Group <- factor(metadata[[sampGroups]])
  title<-paste0('Beta values distribution after ',norm_method,' normalizationn')
  minfi::densityPlot(object, metadata[[sampGroups]],main = title)
  dev.off()
}

# Example
#denplot(filtered,ss,sampGroups,norm_method=norm_function)


