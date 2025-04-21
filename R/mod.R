#' Construct Models and Contrasts with limma
#'
#' @description
#' This function is used to construct models and contrasts for differential expression analysis using limma.
#' It builds a design matrix based on the group variable and optional covariates, fits the linear model, and generates contrasts
#' for pairwise or singular group comparisons.
#'
#' @param object A matrix or data frame containing beta values for the analysis.
#' @param betas_idx A matrix or data frame providing the indices for the beta values if object class is FBM.
#' @param group_var  A string specifying the column name for the grouping variable.
#' @param covs.formula Formula for specifying covariates in the model.
#' @param contrasts Optional string of semicolon-separated contrasts. If not provided, automatic contrasts are generated. A contrast must have the format "Category1-Category2" because the design variable transforms the category values into column names. Then, limma::makeContrasts function interprets the "contrasts" variable and uses it in the model created by the design matrix. 
#' @param covs A character vector specifying additional covariates.
#' @param metadata The metadata or sample sheet.
#' @param set A boolean vector to subset the observations.
#' @param gr A character vector specifying group variables.
#' @param pairwise Logical indicating whether to create pairwise group comparisons (default: TRUE).
#' @param singular Logical indicating whether to create singular group comparisons (default: FALSE).
#' @param rename Character vector to rename the contrasts.
#'
#' @return A limma eBayes model.
#'
#' @export
#'
#' @import limma
#' @import data.table
#' @import bigstatsr
#' @import stringi
#'
#' @examples
#' data("beta_matrix")
#' data("samplesheet")
#' mod(beta_matrix,group_var='Condition',contrasts='Treated-Untreated',metadata=samplesheet)
#'
mod <- function(object, betas_idx = NULL, group_var = "Sample_Group", covs.formula = NULL,
                contrasts = NULL, covs = NULL, metadata, set = TRUE, gr = NULL, pairwise = TRUE,
                singular = FALSE, rename = NULL) {

  # Set row and column names for beta values
  if(class(object)[1] == "FBM"){
    beta_values <- object[]
    rownames(beta_values) <- betas_idx$rn
    colnames(beta_values) <- betas_idx$cn
    object<- beta_values
  }


  # Convert metadata to data.table
  metadata <- data.table::setDT(as.data.frame(metadata))

  # Ensure group_var is a valid variable name
  if(!is.null(group_var)){
    if (!is.numeric(unlist(metadata[,.SD,.SDcols=group_var]))) metadata[, c(group_var) := lapply(.SD, function(x) make.names(x)), .SDcols = c(group_var)]
  }
  # Subset metadata based on 'set'
  metadata <- subset(metadata, set)
  metadata <- droplevels(metadata)

  # Build design matrix
  if (is.null(covs.formula)) {
    if (!is.null(covs) & length(covs) > 0) covs.formula <- paste0("+", paste0(covs, collapse = " + ", sep = ""))

    design <- model.matrix(
      formula(
        paste("~ 0 +", paste0(group_var), covs.formula, sep = " ")
      ),
      data = metadata
    )
  } else {
    design <- model.matrix(formula(covs.formula), data = metadata)
    group_var <- names(attributes(design)$contrasts[1])
  }

  colnames(design) <- make.names(colnames(design))

  # Remove group_var prefix
  out <- tryCatch(
    {
      if (!is.null(rename)) {
        colnames(design) <- rename
      } else {
        colnames(design) <- stringi::stri_replace_all_fixed(colnames(design), group_var, "", vectorize_all = FALSE)
      }
    },
    error = function(cond) {
      colnames(design) <- stringi::stri_replace_all_fixed(colnames(design), group_var, "", vectorize_all = FALSE)
      message(cond)
    }
  )


  # Fit linear model
  fit <- limma::lmFit(object, design)

  # Generate contrasts
  cols <- with(metadata,unique(get(group_var)))
  cont_sing = cont_pair = gr_cont_sing = gr_cont_pair = NULL

  if (pairwise) {
    cont_pair <- apply(combn(cols, 2), 2, function(x) paste(x, collapse = "-"))

    if (!is.null(gr)) {
      gr_cols <- sapply(gr, function(x) contgroup(x, colnames(design)))
      gr_cont_pair <- apply(combn(gr_cols, 2), 2, function(x) paste(x, collapse = "-"))
    }
  }

  if (singular) {
    cont_sing <- apply(combn(cols, length(cols) - 1), 2, function(x) {
      var <- setdiff(cols, x)
      group <- contgroup(group_var, levels = x)
      contrast <- paste0(var, "-", group)
      return(contrast)
    })

    if (!is.null(gr)) {
      gr_cols <- sapply(gr, function(x) contgroup(x, colnames(design)))
      gr_cont_sing <- apply(combn(gr_cols, length(gr_cols) - 1), 2, function(x) {
        var <- setdiff(gr_cols, x)
        group <- contgroup(group_var, levels = x)
        contrast <- paste0(var, "-", group)
        return(contrast)
      })
    }
  }

  cont <- c(cont_sing, cont_pair)
  gr_cont <- c(gr_cont_sing, gr_cont_pair)

  if (is.null(contrasts)) {
    contrasts <- c(cont, gr_cont)
  } else {
    contrasts <- unlist(strsplit(contrasts, ";"))
  }

  contMatrix <- limma::makeContrasts(
    contrasts = contrasts,
    levels = colnames(design)
  )

  # Rename contrasts
  if (!is.null(gr)) colnames(contMatrix) <- stringi::stri_replace_all_fixed(colnames(contMatrix), gr_cols, gr, vectorize_all = FALSE)

  large <- colnames(contMatrix) %in% c(cont_sing, gr_cont_sing)
  colnames(contMatrix)[large] <- sapply(
    colnames(contMatrix)[large], function(x) paste0("sing_", strsplit(x, "-")[[1]][1])
  )

  # Fit model and return eBayes model
  fit2 <- limma::contrasts.fit(fit, contMatrix)
  fit2 <- limma::eBayes(fit2)
  rownames(fit2$cov.coefficients)
  return(fit2)
}

# Helper function to generate group contrasts
contgroup <- function(name, levels) {
  cols <- levels[grepl(name, levels, fixed = TRUE)]
  l <- length(cols)
  paste0("(", paste(cols, collapse = "+"), ")/", l)
}

