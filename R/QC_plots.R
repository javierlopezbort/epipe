#' Generate QC Plots
#'
#' This function generates various quality control (QC) plots for signal distribution to assess the quality of the data before normalization.
#'
#' @param rgSet rgSet object containing channel intenisty values.
#' @param sampGroups A variable to use for coloring groups in the plots. Must match a column name in `rgSet@colData`.
#' @param sampNames A variable containing barcodes or ids matching colnames of the rgsetdata
#' @param pal A color palette to use for group coloring. Defaults to `NULL` (uses the default palette).
#' @param qc_folder A character string specifying the path to the folder where plots will be saved.
#' @param idcol A column name in `rgSet@colData` used as an identifier. Defaults to "barcode".
#'
#'
#' @return @return The function saves QC plots as image files in the specified folder.
#' @details
#' This function generates:
#' - A QC report in PDF format using `minfi::qcReport`.
#' - A density plot showing beta value distributions before normalization.
#' - A mean QC plot showing the mean signal intensities.
#'
#' @export
#' @import minfi
#' @import Biobase
#' @import Biostrings
#' @import S4Vectors
#' @importFrom grDevices png dev.off
#' @importFrom data.table as.data.table
#'
#'
qc<-function(rgSet,sampGroups=NULL, sampNames= "Sample_Name", pal=NULL,
             qc_folder="analysis/intermediate/QC/",idcol="barcode"){

  path=paste0(qc_folder,"Report.pdf")
  colData = rgSet@colData
  ss <- data.table::as.data.table(colData)
  Sample_Group = factor(colData[[sampGroups]])
  if (!(sampNames %in% names(rgSet@colData))) sampNames <- colnames(rgSet)

  # unlist(ss[order(Sample_Group),..idcol])->idx
  # rgSet<-rgSet[,idx]
  # colData = rgSet@colData

  Sample_Group = factor(colData[[sampGroups]])

  # Sample_Group<-setNames(Sample_Group,idx)

  cols<-get_cols(factor_variable = Sample_Group, pal = pal)
  dir.create(qc_folder,recursive=T,showWarnings = F)
  minfi::qcReport(rgSet = rgSet,
                  pdf = paste0(qc_folder,"Report.pdf"),
                  sampGroups= rgSet@colData[[sampGroups]])
                  #sampNames = rgSet@colData[[sampNames]])

  #### Density plot

  grDevices::png(file = paste0(qc_folder, "density_plot_before_norm.png"))
  minfi::densityPlot(rgSet,
                     sampGroups = rgSet@colData[[sampGroups]],main = 'Beta values distribution before normalization')
  grDevices::dev.off()

  # dp <- densPlot(
  #   rgSet, sg = Sample_Group,main = "Beta",
  #   pal =pal#[!duplicated(cols)]
  # )
  # save_plot(dp, path=qc_folder, filename = "density_plot")

  # ggplot2::ggsave(filename = paste0(qc_folder,"density_plot.png"),
  #                 device = "png",plot = dp,
  #                 width = 10, height = 12)


  #### Bean plot

  # grDevices::png(file = paste0(qc_folder,"bean_plot.png"))
  # minfi::densityBeanPlot(
  #   rgSet, sampGroups = rgSet@colData[[sampGroups]],
  #   pal=cols[!duplicated(cols)]
  # )
  # grDevices::dev.off()

  #save_plot(beanplot, path=qc_folder, filename = "bean_plot")


  #### Mean QC plot

  mSet <- minfi::preprocessRaw(rgSet)
  qc   <- minfi::getQC(mSet)
  grDevices::png(file = paste0(qc_folder,"mean_qc.png"))
  minfi::plotQC(qc)
  grDevices::dev.off()
  #save_plot(mean_qc, path=qc_folder, filename = "mean_qc")

}
