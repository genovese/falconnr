#' @useDynLib falconnr
#' @import Rcpp

Rcpp::loadModule("mod_params", TRUE)

# Example use: 
#   p <- LshParameterSetter$new(n, d)
#   p$dump()
#   p$setDistance("negative_inner_product")$setFamily("hyperplane")$dump()


.onLoad <- function(libname, pkgname) {
}

.onUnload <- function (libpath) {
  library.dynam.unload("falconnr", libpath)
}

