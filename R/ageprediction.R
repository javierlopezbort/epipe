#' Age prediction
#'
#' This function predicts biological age from DNA methylation data.
#'
#' @param object MethylSet object
#' @param clock_name Name of the age clock used to predict age. Default=Horvath (horvath 2013). Currently "Hannum", "Levine", "BNN", "skinHorvath", "PedBE", "Wu", "TL", "BLUP", "EN" and "all" are available.
#' @param predict_age Logical value, whether to predict age or not. Default=TRUE
#'
#' @return mSet object with predicted age values in colData().
#' @import methylclock
#' @import minfi
#' @export
#'
#' @examples
#' data("mSet_normalized")
#' ageprediction(mSet_normalized,clock_name='Horvath')

ageprediction<-function(object,clock_name='Horvath',predict_age=TRUE){
  if (predict_age==TRUE){

    beta_values<-minfi::getBeta(object)

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
