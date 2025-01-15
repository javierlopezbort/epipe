#' Summarize DMPs and DMRs
#'
#' This function summarizes Differentially Methylated Positions (DMPs) and Differentially Methylated Regions (DMRs).
#'
#' @param dmps Data.table containing DMP information.
#' @param dmrs Data.table containing DMR information.
#' @param path Path to save the summary, defaults to "results/".
#'
#' @return A data.table summarizing DMPs and DMRs.
#'
#' @import data.table
#' @export

summarize <- function(dmps, dmrs, path = "results/") {
  dir.create(path)
  sdmps <- summary_dmps(dmps, write = FALSE)
  sdmrs <- summary_dmrs(dmrs, write = FALSE)
  s <- merge(sdmps, sdmrs, all = TRUE)
  data.table::fwrite(s, paste0(path, "/summary.csv"))
  return(s)
}
