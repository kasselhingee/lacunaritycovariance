#' @title Evaluate a list of fv objects at a particular argument
#' @export manyfvevalat.R

#' @description Given a list of fv objects evaluates each at the given argument, and returns a list of the same length.
#' @param fvlist A list of fv objects
#' @param argval The value of argument that each function is evaluated at.
#' @param value  The function value to use, ".y" uses the default provided by the fv objects.
#'  Like \code{value} in \code{\link[spatstat]{as.function.fv}} but in this case can only only be a single string.
#' @param extrapolate From \code{\link[spatstat]{as.function.fv}}, whether to extrapolate the fv objects when the provided argval is outside their range.
#' @return A list of values of each fv object at \code{argval}.
#' 

#' @example 
#' fv1 <- Hest(heather$coarse)
#' fv2 <- Hest(complement.owin(heather$coarse))
#' fvlist <- list(fv1,fv2) 

#' fvfunc01 <- as.function.fv(fv1)
#' fvfunc01(0.19)
#' 
#' manyfvevalat(fvlist,0.19)

manyfvevalat <- function(fvlist, argval, value=".y", extrapolate=TRUE){
  stopifnot(length(value)==1)
  fvfuncs <- lapply(fvlist,as.function.fv,value=value,extrapolate=extrapolate)
  fvvals <- mapply(do.call,fvfuncs,args=list(r=list(argval)),SIMPLIFY=FALSE)
  return(fvvals)
}