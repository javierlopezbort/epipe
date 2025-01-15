
#' Apply Noob Normalization to an RGChannelSet Object
#'
#' This function applies the Noob normalization method to an RGChannelSet object.The Noob method (Negative Control
#'  Probe Normalization) is used to preprocess methylation data, adjusting for dye bias in the data.
#'
#' @param rgSet RGChannelSet object to normalize.
#'
#' @return Normalized RGChannelSet object.
#'
#' @importFrom minfi preprocessNoob
#'
#' @export
noob <- function(rgSet){
  minfi::preprocessNoob(rgSet,dyeMethod = 'reference')
}

###############################################################################

#' Apply ssNoob Normalization to an RGChannelSet Object
#'
#' This function applies the ssNoob (Single Sample Negative Control Probe) normalization method to an RGChannelSet object.
#' The ssNoob method is a variant of Noob normalization, designed to address dye bias by using a single sample reference.
#'
#' @param rgSet RGChannelSet object to normalize.
#'
#' @return Normalized RGChannelSet object.
#'
#' @importFrom minfi preprocessNoob
#'
#' @export
ssnoob <- function(rgSet){
   minfi::preprocessNoob(rgSet,dyeMethod = 'single')
}

###############################################################################

#' Apply SWAN Normalization to an RGChannelSet Object
#'
#' This function applies the SWAN (Subset Within Array Normalization) method to an RGChannelSet object.
#' SWAN normalization corrects for dye biases and other systematic errors in methylation data, particularly when working with Illumina's BeadArray platform.
#'
#' @param rgSet RGChannelSet object to normalize.
#'
#' @return Normalized RGChannelSet object.
#'
#' @importFrom minfi preprocessSWAN
#'
#' @export
swan <- function(rgSet){
  minfi::preprocessSWAN(rgSet)
}

###############################################################################

#' Apply FunNorm normalization to RGChannelSet Object
#'
#' This function applies the FunNorm (Functional Normalization) method to an RGChannelSet object. It corrects for biases due to the probe design or other unwanted systematic variations.
#'
#' @param rgSet RGChannelSet object to normalize.
#'
#' @return Normalized RGChannelSet object.
#'
#' @importFrom minfi preprocessFunnorm
#'
#' @export
funn <- function(rgSet){
  minfi::preprocessFunnorm(rgSet)
}

##############################################################################

#' Apply Noob followed by Quantile normalization to RGChannelSet Object
#'
#' @param rgSet RGChannelSet object to normalize.
#'
#' @return Normalized RGChannelSet object.
#'
#' @importFrom minfi preprocessNoob preprocessQuantile
#'
#' @export
noob_pq <- function(rgSet){
  minfi::preprocessNoob(rgSet) |> minfi::preprocessQuantile()
}

##############################################################################

#' Apply Noob followed by SWAN normalization to RGChannelSet Object
#'
#' @param rgSet  RGChannelSet object to normalize.
#'
#' @return Normalized RGChannelSet object.
#'
#' @importFrom minfi preprocessNoob preprocessSWAN
#'
#' @export
noob_swan <- function(rgSet){
  set.seed(123)
  Mset<-minfi::preprocessNoob(rgSet)
  swan<- minfi::preprocessSWAN(rgSet,mSet = Mset)
}


###############################################################################

#' Apply Quantile normalization to RGChannelSetObject
#'
#' @param rgSet RGChannelSet object to normalize.
#'
#' @return Normalized RGChannelSet object.
#'
#' @importFrom minfi preprocessQuantile
#'
#' @export
quantile <- function(rgSet){
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

Em <- function(rgSet, arraytype = NULL){
  require(SummarizedExperiment)
  pd <- rgSet@colData
  qc <- ENmix::QCinfo(rgSet, distplot = FALSE)
  mdat <- ENmix::preprocessENmix(rgSet, QCinfo = qc)
  mSetSqn <- minfi::mapToGenome(mdat)
  mSetSqn@colData <- pd
  return(mSetSqn)
}



#' Normalize Data Using Multiple Methods and Generate Density Plots
#'
#' This function applies various normalization methods to a given dataset and generates density plots
#' for each method. The plots are saved to the specified directory. The user can choose to exclude a specific
#' normalization method by specifying it in the `norm_method` argument.
#'
#'
#' @param object A filtered object (e.g., methylation data) to be normalized
#' @param metadata @param metadata A data frame containing sample metadata, which includes a column for grouping samples (e.g., treatment or control groups).
#' @param sampGroups A character string specifying the column in the metadata that contains the grouping variable for coloring.
#' @param path A character string specifying the path where the density plots will be saved. Default is the current working directory ("./").
#' @param norm_method A character string specifying the normalization method to exclude from the plot generation. Options include 'noob', 'ssnoob', 'swan', 'funn', 'quantile', 'noob_pq', and 'noob_swan'.
#'
#' @import minfi
#' @import grDevices
#'
#' @export
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








