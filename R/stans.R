#' Get a model's full path
#'
#' Code for stan models is stored in the exec/ subdirectory of wherever
#' this package is installed. This is a convenience shortcut to get fully
#' specified paths to a model
#'
#' @param model_name (Relative) filename of model
#' @export
model_filename <- function(model_name) {
  file.path(devtools::inst('beliefupdatr'),
            'exec',
            model_name)
}

#' Compile an included Stan model
#'
#' @param model_name (Relative) filename of model
#' @export
compile_stan <- function(model_name) {
  rstan::stan(model_filename(model_name))
}

