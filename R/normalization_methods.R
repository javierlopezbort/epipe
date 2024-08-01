
#' Apply Noob normalization to RGChannelSet.
#'
#' @param rgSet RGChannelSet object to be normalized.
#' @return Normalized RGChannelSet object.
#'
#' @importFrom minfi preprocessNoob
#' @export
noob <- function(rgSet){
  require(minfi)
  minfi::preprocessNoob(rgSet,dyeMethod = 'reference')
}

###############################################################################

#' Apply ssNoob normalization to RGChannelSet.
#'
#' @param rgSet RGChannelSet object to be normalized.
#' @return Normalized RGChannelSet object.
#'
#' @importFrom minfi preprocessNoob
#' @export
ssnoob <- function(rgSet){
  require(minfi)
  minfi::preprocessNoob(rgSet,dyeMethod = 'single')
}

###############################################################################

#' Apply SWAN normalization to RGChannelSet.
#'
#' @param rgSet RGChannelSet object to be normalized.
#' @return Normalized RGChannelSet object.
#'
#' @importFrom minfi preprocessSWAN
#' @export
swan <- function(rgSet){
  require(minfi)
  minfi::preprocessSWAN(rgSet)
}

###############################################################################

#' Apply FunNorm normalization to RGChannelSet.
#'
#' @param rgSet RGChannelSet object to be normalized.
#' @return Normalized RGChannelSet object.
#'
#' @importFrom minfi preprocessFunnorm
#' @export
funn <- function(rgSet){
  require(minfi)
  minfi::preprocessFunnorm(rgSet)
}

##############################################################################

#' Apply Noob followed by Quantile normalization to RGChannelSet.
#'
#' @param rgSet RGChannelSet object to be normalized.
#' @return Normalized RGChannelSet object.
#'
#' @importFrom minfi preprocessNoob preprocessQuantile
#' @export
noob_pq <- function(rgSet){
  require(minfi)
  minfi::preprocessNoob(rgSet) |> minfi::preprocessQuantile()
}

##############################################################################

#' Apply Noob followed by SWAN normalization to RGChannelSet.
#'
#' @param rgSet RGChannelSet object to be normalized.
#' @return Normalized RGChannelSet object.
#'
#' @importFrom minfi preprocessNoob preprocessSWAN
#' @export
noob_swan <- function(rgSet){
  set.seed(123)
  require(minfi)
  Mset<-minfi::preprocessNoob(rgSet)
  swan<- minfi::preprocessSWAN(rgSet,mSet = Mset)
}


###############################################################################

#' Apply Quantile normalization to RGChannelSet.
#'
#' @param rgSet RGChannelSet object to be normalized.
#' @return Normalized RGChannelSet object.
#'
#' @importFrom minfi preprocessQuantile
#' @export
quantile <- function(rgSet){
  require(minfi)
  minfi::preprocessQuantile(rgSet)
}

##############################################################################

#' Apply ENmix normalization to RGChannelSet.
#'
#' @param rgSet RGChannelSet object to be normalized.
#' @param arraytype Type of array used for preprocessing. If provided, array-specific annotation is used.
#' @return Normalized GenomicRatioSet object.
#'
#' @importFrom ENmix preprocessENmix QCinfo mpreprocess
#' @importFrom minfi makeGenomicRatioSetFromMatrix
#' @export
Em2 <- function(rgSet, arraytype = NULL){
  pd <- rgSet@colData
  qc <- ENmix::QCinfo(rgSet, distplot = FALSE)
  mdat <- ENmix::preprocessENmix(rgSet, QCinfo = qc)
  mSetSqn <- ENmix::mpreprocess(rgSet = rgSet, impute = TRUE)
  if(is.null(arraytype)){
    arraytype <- ifelse(500000 < nrow(mSetSqn), "EPIC", "450K" )
  }

  if(arraytype == "EPIC"){
    mSetSqn <- minfi::makeGenomicRatioSetFromMatrix(
      mat = mSetSqn,
      array = "IlluminaHumanMethylationEPIC",
      annotation = "ilm10b4.hg19"
    )
  } else {
    mSetSqn <- minfi::makeGenomicRatioSetFromMatrix(mSetSqn)
  }

  mSetSqn@colData <- pd
  return(mSetSqn)
}

#' Apply ENmix normalization to RGChannelSet and map to genome.
#'
#' @param rgSet RGChannelSet object to be normalized.
#' @param arraytype Type of array used for preprocessing. If provided, array-specific annotation is used.
#' @return Normalized GenomicRatioSet object.
#'
#' @importFrom ENmix preprocessENmix QCinfo
#' @importFrom minfi mapToGenome
#' @import SummarizedExperiment
#' @export
Em <- function(rgSet, arraytype = NULL){
  require(SummarizedExperiment)
  pd <- rgSet@colData
  qc <- ENmix::QCinfo(rgSet, distplot = FALSE)
  mdat <- ENmix::preprocessENmix(rgSet, QCinfo = qc)
  mSetSqn <- minfi::mapToGenome(mdat)
  mSetSqn@colData <- pd
  return(mSetSqn)
}



#' Function to normalize data by all methods and plot them
#'
#' @param object filtered object to be normalized
#' @param metadata Sample sheet
#' @param sampGroups Sample sheet name variable to color plots
#' @param path  Path to save the plot
#' @param norm_method Normalization method that want to exlcude

normalization_all_methods<-function(object,metadata,sampGroups,path="./",norm_method){

  Sample_Group <- factor(metadata[[sampGroups]])

  # grDevices::pdf(file = paste0(path, "density_plot_all_normalization_methods.pdf"))
  # par(mfrow=c(3,2))

  if (norm_method != 'noob') {
  # Noob
  object.noob<-noob(object)
  grDevices::png(file = paste0(path, "density_plot_noob.png"))
  minfi::densityPlot(object.noob, metadata[[sampGroups]],main = 'Noob',legend = FALSE)
  dev.off()
  }

  if (norm_method != 'ssnoob') {
  # ssnoob
  object.ssnoob<-ssnoob(object)
  grDevices::png(file = paste0(path, "density_plot_ssnoob.png"))
  minfi::densityPlot(object.ssnoob, metadata[[sampGroups]],main = 'ssNoob',legend = FALSE)
  dev.off()
  }


  if (norm_method != 'swan') {
  # Swan
  object.swan<-swan(object)
  grDevices::png(file = paste0(path, "density_plot_swan.png"))
  minfi::densityPlot(object.swan,  metadata[[sampGroups]],main = 'Swan',legend = FALSE)
  dev.off()
  }

  if (norm_method != 'funn') {
  # Funn
  object.funn<-funn(object)
  funn.beta.values<-getBeta(object.funn)
  grDevices::png(file = paste0(path, "density_plot_funn.png"))
  minfi::densityPlot(funn.beta.values,  metadata[[sampGroups]],main = 'Funn',legend = FALSE)
  dev.off()
  }

  if (norm_method != 'quantile') {
  # Quantile
  object.pq<-quantile(object)
  pq.beta.values<-getBeta(object.pq)
  grDevices::png(file = paste0(path, "density_plot_quantile.png"))
  minfi::densityPlot(pq.beta.values,  metadata[[sampGroups]],main = 'Quantile',legend = FALSE)
  dev.off()
  }

  if (norm_method != 'noob_pq') {
  # Noob + quantile
  object.noob.pq<-noob_pq(object)
  noob.pq.beta.values<-getBeta(object.noob.pq)
  grDevices::png(file = paste0(path, "density_plot_noob_quantile.png"))
  minfi::densityPlot(noob.pq.beta.values,  metadata[[sampGroups]],main = 'Noob_Quantile',legend = FALSE)
  dev.off()
  }

  if (norm_method != 'noob_swan') {
  # Noob + swan
  object.noob.swan<-noob_swan(object)
  grDevices::png(file = paste0(path, "density_plot_noob_swan.png"))
  minfi::densityPlot(object.noob.swan, metadata[[sampGroups]],main = 'Noob_swan',legend = FALSE)
  dev.off()
  }
  #dev.off()

}


# Transform norm_function name to a character
makelist<-function(vals){
  vals$norm<-as.character(vals$norm)
  return(vals)

}








