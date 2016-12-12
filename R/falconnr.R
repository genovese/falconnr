#' falconnr: A package for locality-sensitive hashing nearest-neighbor search
#'
#' The falconnr package provides an interface to the
#' \href{https://falconn-lib.org/}{FALCONN} C++ library
#' that supports similarity-search for high-dinensional data
#' based on locality-sensitive hashing.
#' 
#' @author Christopher R. Genovese \email{genovese@@cmu.edu}
#' @references \url{https://falconn-lib.org/}
#'
#' 
#' @docType package
#' @name falconnr
#'
#' @useDynLib falconnr
#' @import Rcpp
#' @importFrom methods new
#' 
NULL

# Module: LSH Contruction Parameter Builder
#   Example use:
#     p <- LshParameterSetter$new(n, d)
#     p$dump()
#     p$distance("negative_inner_product")$family("hyperplane")$dump()

Rcpp::loadModule("mod_params", TRUE)


# Module: LSH Nearest Neighbor Table
#         From FALCONN, this table is static and needs to
#         be reconstructed with each change to the data set.
#
# Example use:
#     nr <- 1000
#     nc <- 10
#     p <- LshParameterSetter$new(nr, nc)
#     X <- matrix(rnorm(nr * nc), nr, nc)
#     tab <- LshNnTable$new(t(X), p)
#     tab$find_nearest_neighbor(as.vector(X[1, ]))

Rcpp::loadModule("mod_table", TRUE)


#' Removes variables loaded from Rcpp modules
#'
#' Because the underlying module pointers exposed by Rcpp
#' do not persist between sessions, is is useful to clean
#' them up before ending the session, lest a later reference
#' cause a crash.
#'
#' This function removes all variables in the calling environment
#' whose class names, as returned by \code{class()}, match
#' the given pattern (\code{"^Rcpp_"} by default).
#'
#' @param  pattern regex matching class of objects to remove.
#'                 Default: match "Rcpp_" at the beginning
#'                 of the class name
#' 
#' @return List of names for removed objects
#' @export
#'
cleanup_modules <- function(pattern="^Rcpp_") {
    objects <- ls(pos=parent.frame())
    class_of <- function(n){class(get(n))}
    removed <- objects[grep(pattern, sapply(objects, class_of))]

    rm(list=removed, pos=parent.frame())
    return( removed )
}               

#' Removes variables and objects loaded from Rcpp modules
#'
#' @seealso \code{\link{falconnr::cleanup_modules}}
#'
#' @export
#' 
cleanup <- function() {
    cleanup_modules("^LshTable$")
    cleanup_modules()
}               


# The following might be unnecessary, but there they are.

.onLoad <- function(libname, pkgname) {
}

.onUnload <- function (libpath) {
  library.dynam.unload("falconnr", libpath)
}

