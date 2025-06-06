#' Apply Filtering to DMRs
#'
#' This function applies filtering to Differentially Methylated Regions (DMRs) and generates plots.
#'
#' @param dmrs Data.table containing DMR information.
#' @param plots Logical, indicating whether to generate plots.
#' @param p.value Vector of p-values.
#' @param mDiff Vector of methylation differences.
#' @param min.cpg Vector of minimum CpG values.
#' @param path Path to save intermediate DMR files.
#'
#' @return A melted data.table containing the filtered DMRs.
#'
#' @import ggplot2
#'
#'

apply_filter_dmrs <- function(dmrs, plots = TRUE, p.value = seq(0.001, 0.11, .025),
                              mDiff = seq(0.15, 0.5, .05), min.cpg = seq(3, 6, 1), path = "analysis/intermediate/dmrs/") {
  if (nrow(dmrs) < 1 | is.null(dmrs)) {

    return(warning("No DMRs available..."))
  } else {
    require(ggplot2)
    dir.create(path)
    params <- expand.grid(p.value, mDiff, min.cpg, s = TRUE)

    res1 <- with(params, Map(function(a, b, c, d) {
      dt <- filter_dmrs(dmrs, a, b, c, d)
      dt$p.val <- a
      dt$mDiff <- b
      dt$min.cpg <- c
      dt
    }, Var1, Var2, Var3, s))
    pdata <- rbindlist(res1)
    colA <- names(pdata)[endsWith(names(pdata), "DMRS")]
    colB <- names(pdata)[endsWith(names(pdata), "Genes")]
    pd <- melt(pdata, measure.vars = list(colA, colB), value.name = c("DMRS", "Genes"),
               variable.name = "type", variable.factor = TRUE)
    levels(pd$type) <- c("Hyper", "Hypo")
    if (plots) {
      plt_list <- list()

      pd2 <- make_ribbon_dt(dt = pd, var = "min.cpg")

      merge(pd, pd2) -> pd3
      (plt_list[["p_Genes"]] <- ggplot2::ggplot(data = pd3, aes(x = mDiff, y = DMRS, group = factor(type), color = factor(type))) +
          geom_line(aes(y = Genes)) +
          theme_minimal() +
          facet_grid(Contrast ~ min.cpg)
      )
      (plt_list[["p_r_Genes"]] <- ggplot2::ggplot(data = pd3, aes(x = mDiff, y = DMRS, group = factor(type), color = factor(type))) +
          geom_line(aes(y = Genes)) +
          geom_ribbon(aes(ymin = min.min.cpg.Genes, ymax = max.min.cpg.Genes,)) +
          facet_grid(Contrast ~ ., margins = FALSE)
      )

      (plt_list[["p_DMRS"]] <- ggplot2::ggplot(data = pd3, aes(x = mDiff, y = DMRS, group = factor(type), color = factor(type))) +
          geom_line() +
          facet_grid(Contrast ~ min.cpg, margins = "min.cpg")
      )
      (plt_list[["p_r_DMRS"]] <- ggplot2::ggplot(data = pd3, aes(x = mDiff, y = DMRS, group = factor(type), color = factor(type))) +
          geom_ribbon(aes(ymin = min.min.cpg.Genes, ymax = max.min.cpg.Genes)) +
          facet_grid(Contrast ~ ., margins = FALSE)
      )

      lapply(1:length(plt_list), function(x)
        save_plot(object = plt_list[[x]],
               filename = paste0(names(plt_list[x]), ".png"),
               path = path
               ))
    }
    return(pd)
  }
}

#' Filter DMRs
#'
#' @title Filter DMRs
#'
#' @description Filters a data table of DMRs based on user-specified thresholds for False Discovery Rate (FDR), methylation difference, and the minimum number of CpGs. Optionally, a summary of the filtered DMRs is generated.
#'
#' @param dmrs A data table containing DMR information.
#' @param p.value Threshold for the False Discovery Rate (default: 0.05).
#' @param mDiff Threshold for the absolute mean methylation difference (default: 0.05).
#' @param min.cpg Minimum number of CpGs in a DMR (default: 3).
#' @param s Logical indicating whether to generate a summary of filtered DMRs (default: FALSE).
#'
#' @return A data table containing the filtered DMRs, or a summary if `s = TRUE`.
#'
#' @import data.table
#'
#' @export
#'

filter_dmrs <- function(dmrs, p.value = 0.05, mDiff = 0.05, min.cpg = 3, s = FALSE) {
  if (length(dmrs) < 1 | is.null(dmrs)) {
    warning("No DMRs available...")
    return(dmrs <- data.table::as.data.table(dmrs))
  } else {

    require(data.table)
    dmrs <- data.table::as.data.table(dmrs)
    #out <- dmrs[abs(meandiff) >= mDiff & no.cpgs >= min.cpg, ]
    out <- dmrs[dmrs$HMFDR <= p.value & abs(dmrs$meandiff) >= mDiff & dmrs$no.cpgs >= min.cpg, ]
    if (s) out <- summary_dmrs(out, write = FALSE)
    return(out)
  }
}

#DMRs
#'
#' @title Summarize DMRs
#'
#' @description Summarizes the information of DMRs, calculating the total number of DMRs, Hyper and Hypo DMRs, and the associated genes for each contrast. Optionally, writes the summary to a specified file.
#'
#' @param dmrs A data table containing DMRs, including information on `Contrast`, `meandiff`, and `overlapping.genes`.
#' @param contrast Optional, a specific contrast to summarize (default: NULL).
#' @param path Path to save the summary file (default: "/results/dmrs/").
#' @param write Logical, whether to save the summary to file (default: TRUE).
#'
#' @return A data table summarizing DMRs with counts of Hyper and Hypo DMRs, and associated genes.
#'
#' @import data.table
#' @import dplyr
#'
#' @export

summary_dmrs <- function(dmrs,contrast=NULL, path = "/results/dmrs/", write = TRUE) {

  if (nrow(dmrs) == 0) {
    # Define the structure of the empty summary dataframe
    summary <- data.table(
      Contrast = unique(dmrs$Contrast),
      Total_DMRS = 0,
      Hyper_DMRS = 0,
      Hypo_DMRS = 0,
      Total_Genes_DMRS = 0,
      Hyper_Genes_DMRS = 0,
      Hypo_Genes_DMRS = 0
    )
    summary <- rbind(summary, list(contrast, 0, 0, 0, 0, 0, 0))

    # Write the empty summary if write is TRUE
    if (write) data.table::fwrite(summary, path)

    return(summary)
  } else {

    dmrs[, Type := ifelse(meandiff > 0, "Hyper", "Hypo")]
    dmrs.l <- dmrs[, .(Total_DMRS=.N,Hyper_DMRS = sum(Type == "Hyper"), Hypo_DMRS = sum(Type == "Hypo")), by = c("Contrast")]
    genes.l <- dmrs[, .(Total_Genes_DMRS=length(unique(unlist(strsplit(overlapping.genes,',')))),Hyper_Genes_DMRS = length(unique(unlist(strsplit(overlapping.genes[Type == "Hyper"], ",")))),
                        Hypo_Genes_DMRS = length(unique(unlist(strsplit(overlapping.genes[Type == "Hypo"], ","))))), by = c("Contrast")]
    summary <- merge(dmrs.l, genes.l)

    # Order the data frame:
    library(dplyr)

    summary <- summary %>%
      select(Contrast, matches("^Total"),matches("Hyper"),matches("Hypo"))


    if (write) data.table::fwrite(summary, path)
  }
  return(summary)
}

# summary_df <- dmrs %>%
#   group_by(Contrast) %>%
#   summarize(
#     Total_DMRs = n(),
#     Total_genes = length(unique(unlist(strsplit(overlapping.genes, ',')))),# Total number of rows per Contrast
#     Hyper_DMRs = sum(Type == 'Hyper'), # Total number of rows with Type 'Hyper'
#     Hyper_genes= length(unique(unlist(strsplit(overlapping.genes[Type == "Hyper"], ",")))),
#     Hypo_DMRs = sum(Type == 'Hypo'), # Total number of rows with Type 'Hypo'
#     Hypo_genes=length(unique(unlist(strsplit(overlapping.genes[Type == "Hypo"], ","))))
#   )
#





#' Make Ribbon Names
#'
#' This function generates ribbon names based on a variable.
#'
#' @param var Variable for which ribbon names are generated.
#'
#' @return A vector of ribbon names.
#'
#' @author Izar de Villasante
#'
#' @examples
#' # Example usage:
#' #ribbon_names <- make_ribbon_names(var = "min.cpg")
#'

make_ribbon_names <- function(var) {
  varnames <- apply(expand.grid(c("min", "max"), var, c("DMRS", "Genes")), 1, function(x) paste(x, collapse = "."))
  return(varnames)
}

#' Make Ribbon Data Table
#'
#' This function generates a data.table for ribbon plotting based on a variable and data.table.
#'
#' @param var Variable for which ribbon data table is generated.
#' @param dt Data.table containing DMR information.
#'
#' @return A data.table for ribbon plotting.
#'
#'
#' @examples
#' # Example usage:
#' #ribbon_dt <- make_ribbon_dt(var = "min.cpg", dt = dt_data)
#'

make_ribbon_dt <- function(var, dt) {
  n <- names(dt)
  sdcols <- setdiff(n, c("Genes", "DMRS", var))
  dt2 <- dt[, make_ribbon_names(var) := .(
    min(DMRS),
    max(DMRS),
    min(Genes),
    max(Genes)
  ), by = sdcols]
  return(dt2)
}
