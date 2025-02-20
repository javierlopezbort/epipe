#' Run EPIPE
#'
#' This function executes the pipeline
#'
#' @import targets
#'
#' @export


run<-function(){
  targets::tar_make()
}

