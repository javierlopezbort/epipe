#' Find Differentially Methylated Regions (DMRs)
#'
#' This function identifies differentially methylated regions (DMRs) from beta values and a linear model.
#'
#' @title Find DMRs
#'
#' @param object DMPextr object or NULL. If NULL, it will be constructed using the model and beta values.
#' @param betas Beta values matrix. Ignored if object is provided.
#' @param model Linear model used for differential methylation analysis.
#' @param fdr False Discovery Rate threshold.
#' @param p.value P-value threshold (not recommended, use fdr).
#' @param bcutoff Beta value cutoff.
#' @param min.cpg Minimum number of CpGs in a DMR.
#' @param ncores Number of CPU cores to use.
#'
#' @return A data.table containing the identified DMRs.
#'
#' @author Izar de Villasante
#' @export
#' @import DMRcate
#' @import S4Vectors
#' @import GenomicRanges
#' @import foreach
#' @import bigstatsr
#' @import data.table
#'
#' @examples
#' # Example usage:
#' model <- lm(betas ~ condition, data = beta_data)
#' dmrs <- find_dmrs(object = NULL, betas = beta_data, model = model)
#'

find_dmrs <- function(object = NULL, betas = NULL, model, fdr = 0.05, p.value = "fdr", bcutoff = 0.3, min.cpg = 5, ncores = NULL) {
  require(DMRcate)
  require(S4Vectors)
  require(GenomicRanges)
  require(foreach)
  require(bigstatsr)
  require(data.table)

  # If beta values are provided and object is not, construct the object using DMPextr
  if (!is.null(betas) & is.null(object)) {
    object <- DMPextr(
      fit = model,
      ContrastsDM = colnames(model$contrasts),
      beta_normalized = betas,
      p.value = 0.95,
      mDiff = 0.01,
      writeOut = FALSE
    )
  }

  conts <- colnames(model$contrasts)

  i = conts

    myAnnotation <- object
    out <- tryCatch(
      {
        object_sub <- object[object$Contrast == i,]

        # Remove chromosomes with 1 or fewer DMPs
        chromosomes <- names(which(table(object_sub$chr) > 1))
        object_sub <- object_sub[object_sub$chr %in% chromosomes,]

        if (nrow(object_sub) < 1) {
          return(data.table(
            seqnames = character(), start = numeric(), end = numeric(), width = numeric(),
            strand = character(), no.cpgs = integer(), min_smoothed_fdr = numeric(),
            Stouffer = numeric(), HMFDR = numeric(), Fisher = numeric(), maxdiff = numeric(),
            meandiff = numeric(), overlapping.genes = character(), Contrast = character()
          ))
        } else {
          annotated <- GenomicRanges::GRanges(
            as.character(object_sub$chr),
            IRanges(object_sub$pos, object_sub$pos),
            stat = object_sub$t,
            diff = object_sub$logFC,
            ind.fdr = object_sub$P.Value,
            is.sig = object_sub$P.Value < 0.05,
            Contrast = i
          )
        }

        names(annotated) <- rownames(object_sub)
        annotated <- sort(annotated)
        myAnnotation <- new("CpGannotated", ranges = annotated)

        DMRs <- DMRcate::dmrcate(
          myAnnotation,
          pcutoff = p.value,  # Use p.value or fdr depending on your preference
          betacutoff = bcutoff,
          min.cpgs = min.cpg,
          C = 2
        )
      },
      error = function(cond) {
        message(cond)
        # Choose a return value in case of error
        new(
          "DMResults",
          coord = character(),
          no.cpgs = integer(),
          min_smoothed_fdr = numeric(),
          Stouffer = numeric(),
          #HMFDR = numeric(),
          Fisher = numeric(),
          maxdiff = numeric(),
          meandiff = numeric()
        ) -> out
      },
      finally = {
        message(paste("Processed Contrast:", i))
        message("Next.")
      }
    )

    if (!identical(out@no.cpgs, integer(0))) {
      results.ranges <- DMRcate::extractRanges(out)
      results.ranges$Contrast = i
    } else {
      return(data.table(
        seqnames = character(), start = numeric(), end = numeric(), width = numeric(),
        strand = character(), no.cpgs = integer(), min_smoothed_fdr = numeric(),
        Stouffer = numeric(), HMFDR = numeric(), Fisher = numeric(), maxdiff = numeric(),
        meandiff = numeric(), overlapping.genes=character(), Contrast = character()
      ))
    }
    data.table::as.data.table(results.ranges)

  results[HMFDR <= fdr, ]  # Filter the results based on the specified FDR threshold
  return(results)
}






