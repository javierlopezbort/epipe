#' Cell type Deconvolution
#'
#' This function performs deconvolution based on DNA methylation data, estimating cell type proportions
#' using raw and normalized methylation datasets. It supports both HM450 and EPIC platforms.
#'
#'
#' @param raw_object The input raw RGChannelSet.
#' @param object MehtylSet object with normalized data.
#' @param arraytype Vector containing the name of the array platform. Options: '450K', 'EPIC', 'EPICv2'.
#'
#' @return A MethylSet object containing deconvolution results in the colData slot.
#'
#' @import FlowSorted.Blood.EPIC
#' @import sesame
#' @import minfi
#'
#' @export
#' @examples
#' data(rgSet)
#' data("mSet_normalized")
#' deconvoluted<-celldeconvolution(rgSet,mSet_normalized,arraytype='EPICv2')
#' colData(deconvoluted)


#Steps:
# 2. Leukocyte proportion
# 3. Cell type proportion estimation


celldeconvolution<-function(raw_object, object, arraytype){

  # Cell type proportion estimation ---> Check issues of immunomethylomics/FlowSorted.Blood.EPIC git repo
  if (arraytype %in% c("EPIC", "450K")){
    cell.type.prop<-estimateCellCounts2(raw_object) # This function includes a preprocess method.
    cell.type.prop<-cell.type.prop$prop

  } else if (arraytype == 'EPICv2'){ # If arraytype is EPICv2 follow a specific workflow
    MSet <-preprocessNoob(raw_object)
    Betas<-getBeta(MSet)
    Betas<- sesame::betasCollapseToPfx(Betas)

    IDOLOptimizedCpGsBloodv2<- FlowSorted.Blood.EPIC::IDOLOptimizedCpGs[which(FlowSorted.Blood.EPIC::IDOLOptimizedCpGs%in%rownames(Betas))]
    cell.type.prop <- projectCellType_CP(
      Betas[IDOLOptimizedCpGsBloodv2, ],
      FlowSorted.Blood.EPIC::IDOLOptimizedCpGs.compTable[IDOLOptimizedCpGsBloodv2,],
      contrastWBC = NULL, nonnegative = TRUE,
      lessThanOne = FALSE
    )

  }

  ## Get beta-values matrix
  beta_values<-minfi::getBeta(object)

  ## Estimate leukocyte fraction (sesame)

  if (arraytype == 'EPICv2'){

    #Be sure that row names are correctly named
    rownames(beta_values)<- sub("_.*$", "", rownames(beta_values))
    leuk.prop<-estimateLeukocyte(beta_values,platform='EPIC')

  } else{
    leuk.prop<-estimateLeukocyte(beta_values,platform=arraytype)
  }


  # Join
  deconvoluted<-cbind(cell.type.prop,leuk.prop)
  dtf<-data.frame(deconvoluted)
  merged_data<-cbind(object@colData,dtf)
  object@colData<-merged_data

  return(object)
}





