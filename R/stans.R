#' @useDynLib beliefupdatr, .registration = TRUE
#' @import Rcpp
#' @import methods
NULL

#' List the included models
#'
#' These can be passed to model_filename or compile_stan
#'
#' @export
list_models <- function() {
  names(stanmodels)
}
