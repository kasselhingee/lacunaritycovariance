#' @title Pointwise summary of a list of im objects
#' @export summary.imlist
#' @method summary imlist
#' @description
#' This function assumes that im objects are each realisations of the same (stochastic) object. 
#' It returns pointwise summaries such as observed sample mean and sample variance.

#' @author Kassel Hingee

#' @param  object A list of im objects
#' @param  harmonizeobject If TRUE (default) the pixel dimensions of the images will be harmonized. Otherwise the object will be tested for compatibility.
#' @param  ... Ignored
#' @return A list im objects containing the pointwise mean, variance and maxima and minima.

#' @examples
#' obspatterns <- replicate(5, rbdd(10, 0.05, window = square(1)), simplify = FALSE)
#' object <- solapply(obspatterns, function(x) balancedracscovariances(x, obswin = square(1), modifications = "pickaH")[[1]])
#' summ <- summary.imlist(object, harmonizeobject = FALSE)

summary.imlist <- function(object, ...,  harmonizeobject = TRUE){
  if (harmonizeobject) {object <- do.call(harmonize.im, args = object)} else { #for pointwise summaries the pixels must represent the same locations for each image.
    stopifnot(do.call(compatible.im, args = object))
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
