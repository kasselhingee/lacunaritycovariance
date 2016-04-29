#' @title Evaluate a list of fv objects at a particular argument
#' @export manyfvevalat.R

#' @description Given a list of fv objects evaluates each at the given argument, and returns a list of the same length.
#' @param fvlist A list of fv objects
#' @return A list of values of each fv object at the given argument.
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
  fvfuncs <- lapply(fvlist,as.function.fv,value=value,extrapolate=extrapolate)
  fvvals <- mapply(do.call,fvfuncs,args=list(r=list(argval)),SIMPLIFY=FALSE)
  
  return(fvvals)
}