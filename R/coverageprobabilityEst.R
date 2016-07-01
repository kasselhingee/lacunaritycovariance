#' @title Estimate the coverage fraction of a stationary RACS
#' 
#' @description Estimates the probability of an aribitray point being inside a stationary RACS.
#' Because the RACS is stationary this is equivalent to estimating \eqn{P(o\in\Xi)}. 
#' The estimate is simply the fraction of \eqn{\Xi} in the observed window.
#' 
#' 
#' @details in reality the area of xi is probably slightly different - depends on the sensing method of xi. Can assume really close though.

#' @param xi An observation of a RACS.
#' @param w The window of observation (not necessarily rectangular)
#' @return an estimate of the coverage probability
#' @author Kassel Hingee 
#' @import spatstat
#' @export coveragefrac
#' @examples
#' xi <- heather$coarse
#' obswindow <- Frame(heather$coarse)
#' coverageProb <- coveragefrac(xi,obswindow)

#' @keywords spatial nonparametric
coveragefrac <- function(xi,w){
  stopifnot(is.owin(xi))
  stopifnot(is.owin(w))   
  xiInsideW <- intersect.owin(xi,w)
  areaxiInside <- area.owin(xiInsideW)
  areaWindow <- area.owin(w)
  
  covProbEstimate <- areaxiInside/areaWindow
  return(covProbEstimate)
}

