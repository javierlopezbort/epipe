#' Get Pathways
#'
#' This function performs pathway enrichment analysis for differentially methylated positions or regions (DMPs or DMRs) using the
#' gprofiler2 package. It processes each contrast defined in the "Contrast" column of the input table and
#' combines the results into a single data.table.
#' It provides pathway analysis for all DMPs, hypermethylated DMPs, and hypomethylated DMPs.
#'
#' @param dmrs A data.table containing Differentially Methylated Positions or Regions (DMPs or DMRs), including at least the 'Gene' and 'Contrast' columns.
#' @param res.folder Path to save the plots and intermediate files.
#' @param cols A character vector specifying the columns to include in the pathway results (default is c("term_size", "query_size", "intersection_size")).
#' @param pval A numeric value for the p-value threshold (default is 0.05).
#' @param topN Number of pathways to show on the summary (default is 50).
#' @param savefile Logical, whether to save the results to a file.
#'
#' @return A data.table containing the pathway results for full, hyper, and hypo DMRs. The table includes columns
#' like 'Contrast', 'FDR', 'p_value', and others specified in the 'cols' parameter. If no data is available,
#' an empty data.table is returned.
#'
#' @import data.table
#' @export
get_pathways <- function(dmrs, res.folder = "results/pathways/", cols = c("term_size", "query_size", "intersection_size"), pval = 0.05, topN = 50, savefile = FALSE) {

  # Convert dmrs to data.table
  data.table::setDT(dmrs)
  varnames = c("Contrast", "FDR", "p_value",cols, "TERM", "method")
  empty_dt <- setNames(data.table(matrix(nrow = 0, ncol = length(varnames))), varnames)
  # Check if there are genes in the input data
  if (nrow(dmrs) < 1) {
    warning("No genes supplied")

    # Create an empty results table with column names

    results <-empty_dt
  } else {
    # Create the results folder if it doesn't exist
    if (savefile == TRUE) suppressWarnings(dir.create(res.folder))

    # Obtain pathway results for all, hyper, and hypo DMRs
    full_pathways <- pathway(dmrs[, .SD, .SDcols = c("Gene", "Contrast")], cols = cols, pval = pval, topN = topN,re.folder=res.folder,enrich_plots=TRUE)
    if(nrow(full_pathways) == 0) full_pathways <- empty_dt
    if (savefile == TRUE) data.table::fwrite(full_pathways, file.path(res.folder, "full_pathway.csv"))

    hyper_pathways <- pathway(dmrs[meandiff > 0, .SD, .SDcols = c("Gene", "Contrast")], cols = cols, pval = pval, topN = topN)
    if(nrow(hyper_pathways) == 0) hyper_pathways <- empty_dt
    if (savefile == TRUE) data.table::fwrite(hyper_pathways, file.path(res.folder, "hyper_pathway.csv"))

    hypo_pathways <- pathway(dmrs[meandiff < 0, .SD, .SDcols = c("Gene", "Contrast")], cols = cols, pval = pval, topN = topN)
    if (savefile == TRUE) data.table::fwrite(hypo_pathways, file.path(res.folder, "hypo_pathway.csv"))
    if(nrow(hypo_pathways) == 0) hypo_pathways <- empty_dt


    # Add a status column to each subset of pathway results
    hypo_pathways$status = "hypo"
    hyper_pathways$status = "hyper"
    full_pathways$status = "both"

    # Combine results for hypo, hyper, and full pathways
    results <- rbind(hypo_pathways, hyper_pathways, full_pathways)
  }

  return(results)
}
