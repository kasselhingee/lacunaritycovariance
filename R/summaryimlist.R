#' @title Pointwise summary of a list of \code{im} objects
#' @export summary.imlist
#' @description
#' This function assumes that \code{im} objects are each realisations of the same (stochastic) object. 
#' It returns pointwise summaries such as observed sample mean and sample variance.

#' @author Kassel Hingee

#' @param  object A list of \code{im} objects
#' @param  harmonizeobject If TRUE (default) the pixel dimensions of the images will be harmonized. Otherwise the object will be tested for compatibility.
#' @param  ... Ignored
#' @return A list \code{im} objects containing the pointwise mean, variance and maxima and minima.

#' @examples
#' # reduce resolution in setcov() for faster (less accurate) computation 
#' oldopt <- spatstat.options()
#' spatstat.options("npixel" = 2^4)
#' 
#' obspatterns <- replicate(3, rbdd(10, 0.05, window = square(1)), simplify = FALSE)
#' ims <- solapply(obspatterns,
#'  function(x) racscovariance(x, obswin = square(1), estimators = "pickaH", drop = TRUE))
#' summ <- summary.imlist(ims, harmonizeobject = FALSE)
#' spatstat.options(oldopt)


#' @export
summary.imlist <- function(object, ...,  harmonizeobject = TRUE){
  if (harmonizeobject) {object <- do.call(harmonize.im, args = object)} else { #for pointwise summaries the pixels must represent the same locations for each image.
    stopifnot(do.call(compatible.im, args = object))
  }
  if ("na.rm" %in% names(list(...)) && list(...)[["na.rm"]]) {
    stop("summary.imlist not able to remove NA values.")
  }
  #class(object) <- "list"
  n <- length(object)
  meanY <- Reduce(Add.im, object)
  meanY <- eval.im(meanY/n)
  sumY2 <- Reduce(Add.im, lapply(object, Square.im))
  varY <- eval.im( (sumY2 - n *( meanY^2))/(n - 1))
  varY <- eval.im(pmax.int(0, varY))
  maxY <- Reduce(Pmax.im, object)
  minY <- Reduce(Pmin.im, object)
  return(solist(mean = meanY, var = varY, max = maxY, min = minY))
}

Square.im <- function(A) { force(A); eval.im(A^2) }
Add.im <- function(A,B){ force(A); force(B); eval.im(A+B, harmonize = FALSE) }

Pmax.im <- function(A, B){force(A); force(B); eval.im(pmax(A, B), harmonize = FALSE)}
Pmin.im <- function(A, B){force(A); force(B); eval.im(pmin(A, B), harmonize = FALSE)}
