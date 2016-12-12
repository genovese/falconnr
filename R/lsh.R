#' A lightweight S4 class representing a FALCONN LSH search table
#'
#' Objects of this class create a Locality-Sensitive Hashing
#' nearest-neighbor search table with an interface to the
#' FALCONN similarity-search library.
#'
#' This is constructed with a data matrix, and the resulting search
#' structure is \emph{static}, in the sense that any change or additino
#' to the data requires that a new LshTable object be created.
#'
#' This is (for now) a lightweight class in that methods for configuring
#' the table and search algorithm are called directly on the slots,
#' which are Rccp objects exposed from C++. See the vignette on
#' configuration for more details. Later versions will expand the range
#' of available methods for this purpose.
#'
#' The main search facilities are encapsulated in a single R method
#' \code{similar}, which can do nearest-neighbor or near-neighbor
#' searches based on the parameters.
#'
#' Note that Rcpp module objects, such as those held in both slots
#' of this class, hold pointers that are not persistent across
#' sessions. To prevent unexpected pointer errors, it is recommended
#' that objects of this class be removed before ending a session.
#' See the function @seealso \code{\link{falconnr::cleanup}} for
#' an automatic way to do this.
#'
#' @slot table  -- a LshNnTable object based on the given data matrix
#' @slot params -- a LshParameterSetter object to configure FALCONN
#' @slot data   -- the data matrix on which the search table was based
#'
#' @export LshTable
#' @exportClass LshTable
LshTable <- setClass("LshTable",
                     slots=c(table="ANY", # Would not recognize "Rcpp_LshNnTable"
                             params="ANY",
                             data="matrix"))

#' Initialize a LshTable given a data matrix
#'
#' @param .Object -- the LshTable object to be initialized
#' @param X       -- the data, a matrix with each \emph{row}
#'                   corresponding to a point.
#'
#' Note the different orientation of \code{X} relative to the
#' interface of the Rcpp-exposed constructor of \code{LshNnTable},
#' which expects each \emph{column} to be a point.
#' 
#' @export
setMethod("initialize",
          signature(.Object="LshTable"),
          function(.Object, X, params=NULL) {
              if ( !is.matrix(X) ) stop("data matrix missing or invalid")
              if ( is.null(params) ) {
                  n <- nrow(X)
                  d <- ncol(X)
                  .Object@params <- LshParameterSetter$new(n, d)
              } else {
                  .Object@params <- params
              }
              
              .Object@table  <- LshNnTable$new(t(X), .Object@params)
              .Object@data   <- X

              return( .Object )
          })

#' Search for data points similar to a given query point
#'
#' @param object -- an object representing data to search
#' @param query  -- a vector whose dimension matches the dimension
#'                  of the data in object
#' 
#' @export
setGeneric("similar", function(object, query, ...) {
                         standardGeneric("similar")
                     })

#' Search for data points close to the given query point
#'
#' @param object -- an LshTable object
#' @param query  -- a vector whose dimension matches the dimension
#'                  of the data in object
#' @param k      -- number of nearest neighbors to find,
#'                  ignored if \code{radius} is supplied
#' @param radius -- a threshold for near-neighbor search
#' @param points -- if FALSE, return indices in the data matrix of found points,
#'                  otherwise, return a submatrix of the found points
#' 
#' @export
#' 
setMethod("similar", "LshTable",
          function(object, query, k=1, radius=NULL, points=FALSE) {
              if ( is.null(radius) ) {
                  if ( k == 1 ) {
                      indices <- object@table$find_nearest_neighbor(query)
                  } else if ( k > 1 ) {
                      indices <- object@table$find_k_nearest_neighbors(query, k)
                  } else {
                      stop("k-nearest-neighbor search for nonpositive k")
                  }
              } else if ( radius < 0.0 ) {
                  stop("near neighbor search with negative radius")
              } else {
                  indices <- object@table$find_near_neighbors(query, radius)
              }

              if ( !points ) {
                  return( indices )
              } else {
                  return( object@X[indices, ] )
              }
          })
