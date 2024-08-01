#' Purify RGChannelSet
#'
#' This function purifies an RGChannelSet by imputing missing values and applying random forest-based correction.
#'
#' @param rgSet RGChannelSet object to be purified.
#' @param arraytype Array type, specifying the type of Illumina array.
#'
#' @return Purified RGChannelSet.
#'
#' @importFrom randomForest randomForest
#' @importFrom impute impute.knn
#' @importFrom cnv.methyl purify
#'
#' @examples
#' # Purify RGChannelSet for EPICv2 array
#' purify(rgSet, arraytype = "EPICv2")
#'
#' # Purify RGChannelSet for EPIC or 450K array
#' purify(rgSet, arraytype = "EPIC")
#'
#' @export
purify <- function(rgSet, arraytype = NULL) {
  library(cnv.methyl)
  if (is.null(arraytype)) {
    guess_arraytype(nrow(rgSet))
  }

  if (arraytype == "EPICv2") {
    library("randomForest")
    library(impute)

    # Load the internal RFpurify_ABSOLUTE function from the cnv.methyl package
    RFpurify_ABSOLUTE <- cnv.methyl:::RFpurify_ABSOLUTE

    # Extract beta values from RGChannelSet
    betas <- minfi::getBeta(rgSet)
    rownames(betas) <- unlist(sapply(strsplit(rownames(betas), "_"), "[")[1, ])

    # Match rows of betas with importance and reorder
    betas <- betas[match(rownames(RFpurify_ABSOLUTE$importance), rownames(betas)), , drop = FALSE]
    rownames(betas) <- rownames(RFpurify_ABSOLUTE$importance)

    # Impute missing values using k-nearest neighbors
    betas <- tryCatch({
      betas <- impute::impute.knn(data = betas, k = 5)$data
    }, error = {
      betas[is.na(betas)] <- mean(betas, na.rm = TRUE)
    })

    # Apply random forest-based correction
    purity <- stats::predict(RFpurify_ABSOLUTE, t(betas))
    colData(rgSet)$purity <- purity

  } else if (arraytype %in% c("EPIC", "450K")) {
    # Use purify function from the cnv.methyl package for EPIC or 450K arrays
    purity <- cnv.methyl::purify(myLoad = rgSet)
    colData(rgSet)$purity <- purity

  } else {
    stop("Unsupported 'arraytype'. Please use 'EPICv2', 'EPIC', or '450K'.")
  }
  return(purity)
}


