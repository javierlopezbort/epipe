#' Generate qc plots
#' @title generate qc plots for signal distribution prior to filtering
#' @param rgSet rgSet object containing channel intenisty values
#' @param sampNames variable containing barcodes or ids matching colnames of the rgsetdata
#' @param sampGroups variables to use for coloring groups
#' @param qc_folder path to the folder where plots will be saved
#' @return plots
#' @author izar de Villasante
#' @export

qc <- function(rgSet,sampGroups=NULL, sampNames= "Sample_Name", cols=NULL,
                qc_folder="analysis/intermediate/QC/"){
  path=paste0(qc_folder,"Report.pdf")
  colData = rgSet@colData
  Sample_Group = colData[[sampGroups]]

  symrgSet = as.symbol(rgSet)
  symsG = as.symbol(Sample_Group)
  symsN = as.symbol(sampNames)
  symqc = as.symbol(qc_folder)
  symcols = as.symbol(cols)
  command_qcReport <- quote(
    minfi::qcReport(rgSet = symrgSet,pdf = path,
                    sampGroups = symsG))
  # sampNames = rgSet@colData[[sampNames]]),
  list(
    tar_target_raw("qc_report",command_qcReport )#,
    # tar_target_raw("KEGG_dict", command_KEGG, deployment = "main"),
    # tar_target_raw("p_dict", quote(rbind(GO_dict,KEGG_dict)),
    #                deployment = "main")
  )
}
