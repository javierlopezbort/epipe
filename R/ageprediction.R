#' Age prediction function
#'
#' This function predicts biological age from DNA methylation data.
#'
#' @param object MethylSet object
#' @param clock_name Name of the age clock used to predict age. Default=Horvath (horvath 2013)
#' @param predict_age Logical value, whether to predict age or not. Default=TRUE
#'
#' @return Data frame with a column containing the predicted age values.
#'
#' @export

# Example:
# clock_name<-'skinHorvath'
# x<-ageprediction(object,clock_name)


# Com a millora es podria posar que es poguessin utilitzar varis clocks alhora.

ageprediction<-function(object,clock_name='Horvath',predict_age=TRUE){
  if (predict_age==TRUE){
    library(methylclock)
    beta_values<-getBeta(object)

    #Be sure that row names are correctly named (specially for EPICv2 arrays)
    rownames(beta_values)<- sub("_.*$", "", rownames(beta_values))

    #Predict age
    methylclock<-DNAmAge(beta_values,clocks = clock_name)
    predicted_age<-methylclock[,2]

    #Add the predicted age to the object metadata
    object$predictedAge<-as.numeric(unlist(predicted_age))

  }

  return(object)
}

#clean<-ageprediction(clean_EPICv2)

