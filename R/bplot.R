#' Generate PCA bi-plots
#'
#' This function generates PCA bi-plots using ggplot2 and related packages.
#'
#' @title Generate PCA bi-plots
#'
#' @param pca A prcomp object (output of PCA analysis).
#' @param ss Data frame or data table containing samplesheet with metadata information.
#' @param colgroup Character column name in ss. Color points based on this variable.
#' @param s Character column name in ss. Shape points according to categories of that variable.
#' @param combs 2D matrix (default = combn(4,2)) mapping PCs to plot (1st row = X, 2nd row = Y).
#' @param cols Color palette to use (default: get_cols() function).
#' @param tit Title for the plot.
#' @param labs Logical. Add labels if TRUE (default).
#' @param overlap Numeric (default = Inf). Pass to max.overlaps in geom_text_repel.
#' @param alfa Alpha value for transparency (default = 0.3).
#' @param folder Path to folder where plots are saved (default: "analysis/pca/bplots/").
#' @param idcol Column name for the sample ID (default: "Sample_Name").
#'
#' @return Path to the folder where results are saved.
#'
#' @author izar de Villasante
#' @export
#'
#' @import ggplot2
#' @import gplots
#' @import ggrepel
#' @import ggfortify
#' @import RColorBrewer
#'
bplot <- function(pca, ss, colgroup, s, combs = NULL, pal=NULL, tit = NULL, labs = TRUE, overlap = Inf, alfa = 0.3, folder = "analysis/pca/bplots/", idcol = "Sample_Name") {

  # Check for required packages and install if necessary
  # install_if_missing <- function(pkg) {
  #   if (!requireNamespace(pkg, quietly = TRUE)) renv::install(pkg)
  # }
  # sapply(c("ggplot2", "gplots", "ggrepel", "ggfortify", "RColorBrewer"), install_if_missing)

  require(ggplot2)
  require(gplots)
  require(ggrepel)
  require(ggfortify)
  require(RColorBrewer)

  # Ensure ss is a data.table
  ss <- droplevels.data.frame(ss)
  ss <- data.table::as.data.table(ss)
  ss <- ss[get(idcol) == rownames(pca$x)]
  n <- nrow(pca$rotation)

  # Iterate over colgroup variables
  lapply(colgroup, function(f) {

    # Set color scale based on whether the variable is numeric or categorical
    if (is.numeric(ss[[f]])) {
      myPalette <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))
      sc <- scale_colour_gradientn(colours = myPalette(100), limits = range(ss[[f]]))
    } else {
      cols<-get_cols(factor_variable = factor(ss[[f]]), pal = pal)
      sc <- scale_color_manual(values = unique(cols))
    }

    # Set default combs if NULL
    if (is.null(combs)) combs <- combn(4, 2)

    # Iterate over combs to create plots
    for (i in 1:dim(combs)[2]) {
      if (is.null(tit)) tit <- paste0("colored.by.", f, "_shape.", s)
      ap <- ggplot2::autoplot(pca, x = combs[1, i], y = combs[2, i], data = ss, colour = f, shape = s, alpha = alfa, size = 1) +
        sc +
        ggtitle(tit) +
        theme_bw(base_size = 7) +
        theme(legend.key = element_blank(), legend.key.size = unit(1, "point"))

      # Add labels if labs is TRUE
      if (isTRUE(labs)) {
        ap <- ap +
          geom_text_repel(
            aes(label = Sample_Name, color = with(ss, get(f))),
            show.legend = FALSE, size = 1.5, max.overlaps = overlap, segment.size = 0.2,
            min.segment.length = 0.8, point.size = 0.5
          ) +
          labs(colour = f)
      }

      # Create folder if not exists and save plot
      dir.create(folder)
      ggsave(paste0(folder, tit, i, "_", n, ".png"), plot = ap, width = 5.56, height = 2.80, units = "in")
    }
  })

  return(folder)
}
