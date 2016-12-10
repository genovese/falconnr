#' @useDynLib falconnr
#' @import Rcpp

# Module: LSH Contruction Parameter Builder
#   Example use:
#     p <- LshParameterSetter$new(n, d)
#     p$dump()
#     p$distance("negative_inner_product")$family("hyperplane")$dump()

Rcpp::loadModule("mod_params", TRUE)


# The following might be unnecessary, but there they are.

.onLoad <- function(libname, pkgname) {
}

.onUnload <- function (libpath) {
  library.dynam.unload("falconnr", libpath)
}

