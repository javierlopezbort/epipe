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
#' @param pal Color palette to use (default: get_cols() function).
#' @param tit Title for the plot.
#' @param labs Logical. Add labels if TRUE (default).
#' @param overlap Numeric (default = Inf). Pass to max.overlaps in geom_text_repel.
#' @param alfa Alpha value for transparency (default = 0.3).
#' @param folder Path to folder where plots are saved (default: "analysis/pca/bplots/").
#' @param idcol Column name for the sample ID (default: "Sample_Name").
#'
#' @return Path to the folder where results are saved.
#'
#' @export
#'
#' @import ggplot2
#' @import gplots
#' @import ggrepel
#' @import ggfortify
#' @import RColorBrewer
#'
#' @examples
#' \dontrun{
#' data("top_500_beta")
#' data("samplesheet")
#' pca<-pca_res(top_500_beta,samplesheet,sampGroups='Condition',scale = TRUE, center = TRUE,filename='pca_plot')
#'
#' bplots_var<-c('Sentrix_ID','Sentrix_Position','Condition')
#' bplot(pca,samplesheet,colgroup=bplots_var,s='Condition',folder='prova/')
#' }

bplot <- function(pca, ss, colgroup, s, combs = NULL, pal = NULL, tit = NULL, labs = TRUE, overlap = Inf, alfa = 0.3, folder = "analysis/pca/bplots/", idcol = "Sample_Name") {

  # Ensure ss is a data.table
  ss <- droplevels.data.frame(ss)
  ss <- data.table::as.data.table(ss)
  ss <- ss[get(idcol) == rownames(pca$x)]
  n <- nrow(pca$rotation)

  # Iterate over colgroup variables
  lapply(colgroup, function(f) {

    # Set color scale based on whether the variable is numeric or categorical
    if (is.numeric(ss[[f]])) {
      # Define custom color palette excluding middle light colors
      customPalette <- function(n) {
        allColors <- rev(RColorBrewer::brewer.pal(11, "RdYlGn"))
        # Exclude the middle light colors (e.g., the middle 3 colors)
        limitedColors <- allColors[c(1:4, 8:11)]
        colorRampPalette(limitedColors)(n)
      }
      sc <- scale_colour_gradientn(
        colours = customPalette(100),
        limits = range(ss[[f]]),
        guide = guide_colorbar(
          barwidth = 0.8,
          barheight = 8,
          title.position = "top",
          title.hjust = 0.5,
          label.position = "right",
          ticks.colour = "black",
          ticks.linewidth = 3
        )
      )
    } else {
      cols <- get_cols(factor_variable = factor(ss[[f]]), pal = pal)
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
        theme(
          legend.key = element_blank(),
          legend.key.size = unit(1, "point"),
          legend.text = element_text(size = 5),  # Increase legend text size
          legend.title = element_text(size = 5),  # Increase legend title size
          legend.position = "right",  # Position the legend on the right
          legend.margin = margin(6, 6, 6, 6)
        )

      # Add labels if labs is TRUE
      if (isTRUE(labs)) {
        ap <- ap +
          geom_text_repel(
            aes(label = get(idcol), color = with(ss, get(f))),
            show.legend = FALSE, size = 1.5, max.overlaps = overlap, segment.size = 0.2,
            min.segment.length = 0.8, point.size = 0.5
          ) +
          labs(colour = f)
      }

      # Create folder if not exists and save plot
      dir.create(folder, showWarnings = FALSE, recursive = TRUE)
      ggsave(paste0(folder, tit, i, "_", n, ".png"), plot = ap, width = 5.56, height = 2.80, units = "in")
    }
  })

  return(folder)
}
