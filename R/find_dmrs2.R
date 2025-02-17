#' Find Differentially Methylated Regions (DMRs) -2
#'
#' This function identifies differentially methylated regions (DMRs) from beta values and a linear model.
#'
#' @title Find DMRs
#'
#' @param beta_values Beta values matrix.
#' @param model Linear model used for differential methylation analysis.
#' @param fdr False Discovery Rate threshold.
#' @param p.value P-value threshold (not recommended, use fdr).
#' @param bcutoff Beta value cutoff.
#' @param min.cpg Minimum number of CpGs in a DMR.
#' @param output_dir Path to folder to save results.
#' @param arraytype Type of array. EPICv2, EPIC or 450K.
#
#'
#' @return A data.table containing the identified DMRs.
#'
#' @export
#' @import DMRcate
#' @import S4Vectors
#' @import GenomicRanges
#' @import data.table
#'
#'
# @examples
# # Example usage:
# data(beta_matrix)
# data(model)
# dmrs<-find_dmrs(beta_values = beta_matrix,model = model,pal=c("#1B9E77","#7570B3","#D95F02"),output_dir = './',arraytype='EPICv2')


find_dmrs<-function(beta_values,model,fdr = 0.01, p.value = "fdr", bcutoff = 0.05, min.cpg = 5,pal,output_dir='./',arraytype){

  m_values <- minfi::logit2(beta_values)

  # Contrast matrix
  contrast_matrix <- model$contrasts

  # Create colors for the plots
  sample_type <- apply(model$design, 1, function(row) {
    if (row["Control"] == 1) {
      return("Control")
    } else if (row["Treated"] == 1) {
      return("Treated")
    } else if (row["Untreat"] == 1) {
      return("Untreat")
    }
  })


  sample_conditions <- colnames(model$design)
  group_colors<-setNames(pal,sample_conditions)
  cols <- group_colors[sample_type]

  final_dmr_table<-data.table::data.table()

  # Find DMRs for each contrast
  for (contrast in colnames(contrast_matrix)) {

    if (arraytype=='EPICv2'){
      my_annotation <- cpg.annotate(
        datatype = 'array',
        object = m_values,
        what = 'M',
        analysis.type = "differential",
        arraytype = "EPICv2",
        design = model$design,
        contrasts = TRUE,
        cont.matrix = contrast_matrix,
        epicv2Filter = "mean",fdr=fdr,
        coef = contrast)
    } else{
      my_annotation <- cpg.annotate(
        datatype = 'array',
        object = m_values,
        what = 'M',
        analysis.type = "differential",
        arraytype = arraytype,
        design = model$design,
        contrasts = TRUE,
        cont.matrix = contrast_matrix,
        fdr=fdr,
        coef = contrast)

    }


    message(paste('contrast:',contrast))

    element_metadata <- my_annotation@ranges@elementMetadata

    if (any(element_metadata$is.sig == TRUE)) { # if there are TRUEs (significant CPG positions)
      dmrcoutput<-dmrcate(my_annotation,
                          lambda = 1000,
                          # pcutoff = p.value,
                          betacutoff = bcutoff,
                          min.cpgs = min.cpg,
                          C = 2)


      # Check if there are DMRs before extracting ranges
      if (length(dmrcoutput@coord) == 0) {
        dmr_dt <- data.table::data.table(
          seqnames = character(), start = numeric(), end = numeric(), width = numeric(),
          strand = character(), no.cpgs = integer(), min_smoothed_fdr = numeric(),
          Stouffer = numeric(), HMFDR = numeric(), Fisher = numeric(), maxdiff = numeric(),
          meandiff = numeric(), overlapping.genes = character(), Contrast = character()) # EMPTY DT

      } else{
        results.ranges<-extractRanges(dmrcoutput,genome='hg38')

        dmr_dt<-data.table::data.table(as.data.frame(results.ranges))
        dmr_dt<-dmr_dt[,Contrast:=contrast]

        # Sort DMRS
        sorted_results <- results.ranges[order(results.ranges$min_smoothed_fdr, -abs(results.ranges$meandiff))]

        # Select top 10 DMRS for the plot (if there are up to 10)
        num_top_dmrs<-min(10, length(sorted_results))
        top_dmrs <- sorted_results[1:num_top_dmrs]

        # for (i in 1:num_top_dmrs) {
        #
        #   plot_file <- paste0(output_dir, contrast, "_", i, ".png")
        #   png(plot_file, width = 800, height = 600)
        #
        #   DMR.plot(ranges = results.ranges,dmr = i,CpGs = my_annotation, what = "Beta",
        #            arraytype = "EPICv2",phen.col = cols,genome = "hg38")
        #
        #   dev.off()
        #
        # }
      }


    }

    else{
      dmr_dt <- data.table::data.table(
        seqnames = character(), start = numeric(), end = numeric(), width = numeric(),
        strand = character(), no.cpgs = integer(), min_smoothed_fdr = numeric(),
        Stouffer = numeric(), HMFDR = numeric(), Fisher = numeric(), maxdiff = numeric(),
        meandiff = numeric(), overlapping.genes = character(), Contrast = character()) # EMPTY DT
    }

    final_dmr_table<-rbind(final_dmr_table,dmr_dt,fill=TRUE)

  }

  return(final_dmr_table)
}
