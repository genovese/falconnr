#' @useDynLib falconnr
#' @import Rcpp

# Example use of mod_params: 
#   p <- LshParameterSetter$new(n, d)
#   p$dump()
#   p$distance("negative_inner_product")$family("hyperplane")$numHashFunctions(3)$numHashTables(2)$dump()

Rcpp::loadModule("mod_params", TRUE)


.onLoad <- function(libname, pkgname) {
}

.onUnload <- function (libpath) {
  library.dynam.unload("falconnr", libpath)
}

