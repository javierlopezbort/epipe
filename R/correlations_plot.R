#' Correlation analysis plot
#'
#' @description Pairwise correlation plots. Plots are colored if the variables are correlated.
#'
#' @param correlation_list A list of correlated variable pairs, where each entry is a vector containing:
#'   - The first variable name.
#'   - The second variable name.
#'   - The correlation score or p-value as a string (e.g., "0.05").
#'
#' @param dataset A data frame containing the sample sheet with variables to analyze.
#' @param sampGroups A string specifying the grouping variable for coloring points.
#' @param path A string specifying the path to save the correlation plot (default: './').
#'
#' @import GGally
#' @import ggplot2
#' @import rlang
#'
#' @return Saves the correlation plot as a PNG file in the specified path.
#'


# List with variable correlations
#correlation_list<-list(c('Sample_Group','Sentrix_ID'),c('Condition','Sentrix_ID'),c('Type','arraytype'))


correlation_plot<-function(correlation_list,dataset,sampGroups,path = "./"){

  # Initialize ggpairs plot
  pm <- ggpairs(dataset, upper = 'blank')

  for (i in seq_along(correlation_list)) {

    positions<-which(colnames(dataset) %in% correlation_list[[i]])


    # Extract plot coordinates for current correlation pair
    p <- pm[positions[[2]], positions[[1]]]

    #specify variables
    var1<-correlation_list[[i]][[1]]
    var2<-correlation_list[[i]][[2]]

    is_var1_categorical <- is.factor(dataset[[var1]])
    is_var2_categorical <- is.factor(dataset[[var2]])

    if (is_var1_categorical & !is_var2_categorical) {
      select_variable<-var1

    } else if (!is_var1_categorical & is_var2_categorical) {
      select_variable<-var2

    } else if (is_var1_categorical & is_var2_categorical){
      select_variable<-var2

    } else {
      select_variable<-'red'
    }

    if (select_variable=='red'){
      #print(select_variable)
      p <- p + geom_point(aes(color = !!rlang::sym(sampGroups)))
      p<-p+ggplot2::theme(plot.background = ggplot2::element_rect(colour = 'red',size=2))

      # Get correlation score
      correlation_score<-round(as.numeric(correlation_list[[i]][[3]]),3)
      t<-ggally_text(
        label=paste('Correlation!\n','Cor.Score:',correlation_score),
        col='red',size=3)+theme_classic()

    }else{
      # Add aesthetics
      p <- p + aes(fill = !!rlang::sym(sampGroups))
      #p <- p + aes(fill = !!rlang::sym(select_variable))###this variable can change or can be specified by a color ex: 'red'
      p<-p+ggplot2::theme(plot.background = ggplot2::element_rect(fill = 'red',size=3))

      # If variables are categorical-->get the p.value of the test and add it to the graphic.
      p.value<-round(as.numeric(correlation_list[[i]][[3]]),3)

      t<-ggally_text(
        label=paste('Correlation!\n','P.value:',p.value),
        col='red',size=3)+theme_classic()
    }

    # Update plot for coloring correlations
    pm[positions[[2]], positions[[1]]] <- p

    # Update plot for text mark
    pm[positions[[1]], positions[[2]]] <- t

  }


  # Save plot:

  png(file=file.path(path,paste0('correlation_analysis.png')),width=1000,height = 1000)
  print(pm)
  dev.off()

  #save_plot(pm,filename = 'Correlation_analysis', path=path)
  #return (pm)
}





