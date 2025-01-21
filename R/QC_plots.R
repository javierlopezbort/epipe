#' Generate QC Plots
#'
#' This function generates various quality control (QC) plots for signal distribution to assess the quality of the data before normalization.
#'
#' @param rgSet rgSet object containing channel intenisty values.
#' @param sampGroups A variable to use for coloring groups in the plots. Must match a column name in `rgSet@colData`.
#' @param idcol A column name in `rgSet@colData` used as an identifier. Defaults to "Sample_Name".
#' @param pal A color palette to use for group coloring. Defaults to `NULL` (uses the default palette).
#' @param qc_folder A character string specifying the path to the folder where plots will be saved.
#'
#'
#'
#' @return The function saves QC plots as image files in the specified folder.
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
#'
#' @examples
#'
#' data("rgSet")
#'
#' qc(rgSet,
#' sampGroups = 'Type',
#' idcol='Sample_Name',
#' qc_folder = './')
#'
#'
qc<-function(rgSet,sampGroups=NULL, idcol="Sample_Name", pal=NULL,
             qc_folder="analysis/intermediate/QC/"){

  colData = rgSet@colData
  Sample_Group = factor(colData[[sampGroups]])

  cols<-get_cols(factor_variable = Sample_Group, pal = pal)
  dir.create(qc_folder,recursive=T,showWarnings = F)

  minfi::qcReport(rgSet = rgSet,
                  pdf = paste0(qc_folder,"Report.pdf"),
                  sampGroups= Sample_Group)


  #### Density plot

  grDevices::png(file = paste0(qc_folder, "density_plot_before_norm.png"))
  minfi::densityPlot(rgSet,
                     sampGroups =Sample_Group,main = 'Beta values distribution before normalization')
  grDevices::dev.off()


  #### Mean QC plot

  mSet <- minfi::preprocessRaw(rgSet)
  qc   <- minfi::getQC(mSet)
  grDevices::png(file = paste0(qc_folder,"mean_qc.png"))
  minfi::plotQC(qc)
  grDevices::dev.off()


}
