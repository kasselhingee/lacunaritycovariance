#' @title Pointwise summary of a list of im objects
#' @export summary.imlist
#' @description
#' This function assumes that im objects are each realisations of the same (stochastic) object. 
#' It returns pointwise summaries such as observed sample mean and sample variance.

#' @author Kassel Hingee

#' @param  ims A list of im objects
#' @param  harmonizeims If TRUE (default) the pixel dimensions of the images will be harmonized. Otherwise the ims will be tested for compatibility.
#' @return A list im objects containing the pointwise mean, variance and maxima and minima.

#' @examples
#' obspatterns <- replicate(5, rbdd(10, 0.05, window = square(1)), simplify = FALSE)
#' ims <- solapply(obspatterns, function(x) balancedracscovariances(x, obswin = square(1), modifications = "pickaH")[[1]])
#' summ <- summary.imlist(ims, harmonizeims = FALSE)

summary.imlist <- function(ims, harmonizeims = TRUE){
  if (harmonizeims) {ims <- do.call(harmonize.im, args = ims)} else { #for pointwise summaries the pixels must represent the same locations for each image.
    stopifnot(do.call(compatible.im, args = ims))
  }
  #class(ims) <- "list"
  n <- length(ims)
  meanY <- Reduce(Add, ims)
  meanY <- eval.im(meanY/n)
  sumY2 <- Reduce(Add, lapply(ims, Square))
  varY <- eval.im( (sumY2 - n *( meanY^2))/(n - 1))
  varY <- eval.im(pmax.int(0, varY))
  maxY <- Reduce(Pmax, ims)
  minY <- Reduce(Pmin, ims)
  return(solist(mean = meanY, var = varY, max = maxY, min = minY))
}

Square <- function(A) { force(A); eval.im(A^2) }
Add <- function(A,B){ force(A); force(B); eval.im(A+B, harmonize = FALSE) }

Pmax <- function(A, B){force(A); force(B); eval.im(pmax(A, B), harmonize = FALSE)}
Pmin <- function(A, B){force(A); force(B); eval.im(pmin(A, B), harmonize = FALSE)}
