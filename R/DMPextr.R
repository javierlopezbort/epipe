#' For each contrast extract a set of DMPs and add gene annotation and methylation values
#'
#'
#' @title Extract DMPs, annotation and methylation difference for each contrast
#'
#' @return data.table
#' @author Izar de Villasante
#' @export
#' @import minfi
#' @import data.table
#' @import limma
#' @param beta_normalized normalized betavalues, as produce by minfi::getBeta(grSet_noob)),
#'  where colnames(beta_normalized) == metadata$sample_Name
#' @param ContrastsDM list of contrasts as returned by limma::makeContrasts()
#' which will pass to limma topTable as input
#' @param mDiff absolute mean methylation difference between groups to filter by
#' @param ann annotation dataset from manifest with metadata such as gene info,
#' CGI, RefGene, etc. see topTable genelist arg.
#' @param writeOut save result as .csv default = TRUE.
#' @param writedir
#'
#' @inheritParams limma::topTable
#' @examples
#'
#' betas<-readRDS("data/beta_noob.rds")
#' fit<-readRDS("data/fit2.rds")
#' ann<-readRDS("data/ann.rds")
#' DMPann <- DMPextr(fit = fit,                       # linear contrast model
#'                   ContrastsDM = ContrastsDM,          # contrasts
#'                   p.value = 0.01,                      # filter significantly different probes
#'                   beta_normalized = beta_noob,        # extract mean group betas
#'                   mDiff = 0.5,                        # select mean methylation differences
#'                   ann = ann,                          # annotate positions (CGI, RefGene, etc)
#'                   writeOut = FALSE                    # write output to file
DMPextr <- function(
    fit, ContrastsDM = colnames(fit$contrasts), p.value, beta_normalized, betas_idx = NULL,
    mDiff, ann = NULL, writedir = "analysis/DMP_", writeOut = TRUE, ncores = NULL, columns = TRUE
) {
  require(foreach)
  require(bigstatsr)
  require(data.table)

  # Check if annotation is provided, if not, use Illumina annotation
  if (is.null(ann)) {
    dt_epic <- as.data.table(IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Other, keep.rownames = "ProbeID")
    dt_450k <- as.data.table(IlluminaHumanMethylation450kanno.ilmn12.hg19::Other, keep.rownames = "ProbeID")
    ann <- merge(dt_epic, dt_450k, all = TRUE)
  }


  # Create a bigstatsr FBM (File-backed Big Matrix) for beta values
  if(class(beta_normalized)[1] != "FBM"){
    ann <- ann[rownames(beta_normalized), ]
    betas <- bigstatsr::FBM(NROW(beta_normalized), NCOL(beta_normalized))
    betas[] <- as.matrix(beta_normalized)
  }else{
    betas <- beta_normalized
    ann <- ann[betas_idx$rn, ]
  }

  ann <- data.table::as.data.table(ann, keep.rownames = "ProbeID")
  data.table::setkey(ann, "ProbeID")

  ann$rn <- 1:nrow(ann)

  # Set the number of cores for parallel processing
  if (is.null(ncores)) ncores <- length(ContrastsDM)

  # Initialize parallel cluster
  cl <- parallel::makeCluster(ncores, outfile = "", useXDR = FALSE, type = "FORK")
  parallel::clusterEvalQ(
    cl, {
      requireNamespace(c("limma", "data.table"))
      data.table::setDTthreads(0)
    }
  )
  doParallel::registerDoParallel(cl)

  message("Processing ", length(ContrastsDM), " contrasts. Using ", ncores, " cores.")

  # Use parallel processing to iterate through contrasts
  res <- foreach::foreach(
    i = itertools::isplitIndices(length(ContrastsDM), chunks = ncores),
    .combine = 'rbind',
    .packages = c("limma", "data.table"),
    .inorder = FALSE,
    .errorhandling = "pass"
  ) %dopar% {
    # Extract DMPs using limma::topTable
    DMP_1 <- limma::topTable(
      fit,
      num = Inf,
      coef = i,
      genelist = ann,
      p.value = p.value
    )
    DMP_1$rn<-1:nrow(DMP_1)
    # If no DMPs found, return a data.table with proper column names
    if (nrow(DMP_1) < 2) {
      warning(paste("No DMP found for contrast:", ContrastsDM[i], sep = " "))
      dt <- ann[0, ]
      dt[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B") := numeric()]
      dt$Type <- character(length = 0L)
      dt$Contrast <- character(length = 0L)
      dt$diff_meanMeth <- numeric(length = 0L)
      return(dt)
    } else {
      # Add a "Type" column based on logFC
      DMP_1$Type <- "Hyper"
      DMP_1$Type[which(DMP_1$logFC < 0)] <- "Hypo"

      # Extract specific columns from annotation
      # DMP_1[, columns]

      # Add "Contrast" column
      DMP_1$Contrast <- ContrastsDM[i]

      # Extract methylation values for DMPs and samples
      c1 <- fit$contrasts[, ContrastsDM[i]]
      c1 <- c1[c1 != 0]
      vars1 <- names(c1)
      design <- fit$design[, vars1]
      design <- t(design) * c1

      b <- betas[DMP_1$rn, ]
      DMP_1$diff_meanMeth <- rowSums(
        apply(
          design,
          1,
          function(x) apply(b %*% diag(x), 1, function(y) mean(y[y != 0]))
        )
      )

      # Filter for absolute methylation difference
      DMP_1 <- DMP_1[which(abs(DMP_1$diff_meanMeth) >= mDiff), ]

      # Write output file
      if (writeOut == TRUE) {
        cat(paste("writing analysis/DMP_", ContrastsDM[i], ".csv\n", sep = ""))
        data.table::fwrite(DMP_1, file = paste(writedir, ContrastsDM[i], ".csv", sep = ""))
      }

      return(DMP_1)
    }
  }

  # Stop parallel cluster
  parallel::stopCluster(cl)
  return(res)
}
