#' @title Estimate the coverage fraction of a stationary RACS
#' 
#' @description 
#' The coverage fraction of a stationary RACS is also known as the coverage probability and is the probability that an arbitrary point is covered by the RACS.
#' Given a realisation, \eqn{X}, in window \eqn{W} of a stationary RACS the coverage fraction estimate that this function calculates is the fraction of \eqn{W} covered by \eqn{X}.
#' This estimate relies on the window \eqn{W} being large compared to the spatial dependence of the RACS.
#' See [1, section 6.4.2] for more theoretical details.
#' 
#' @references [1] Chiu, S.N., Stoyan, D., Kendall, W.S. and Mecke, J. (2013) Stochastic Geometry and Its Applications, 3rd ed. Chichester, United Kingdom: John Wiley & Sons.
 
#' 
#' @param xi An observation of a RACS in \pkg{spatstat}'s \code{owin} format.
#' @param w The window of observation (not necessarily rectangular) also in \code{owin} format.
#' @return An estimate of the coverage probability
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

