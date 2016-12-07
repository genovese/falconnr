#' @useDynLib falconnr
#' @importFrom Rcpp sourceCpp


.onUnload <- function (libpath) {
  library.dynam.unload("falconnr", libpath)
}

