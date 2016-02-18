#' @title Estimate the coverage fraction of a stationary RACS
#' 
#' @description Estimates the probability of an aribitray point being inside a stationary RACS.
#' Because the RACS is stationary this is equivalent to estimating \eqn{P(o\in\Xi)}. 
#' The estimate is simply the fraction of \eqn{\Xi} in the observed window.
#' 
#' 
#' @details in reality the area of Xi is probably slightly different - depends on the sensing method of Xi. Can assume really close though.

#' @param Xi An observation of a RACS.
#' @param w The window of observation (not necessarily rectangular)
#' @return an estimate of the coverage probability
#' @author Kassel Hingee 
#' @import spatstat
#' @export coveragefrac
coveragefrac <- function(Xi,w){
  stopifnot(is.owin(Xi))
  stopifnot(is.owin(w))   
  XiInsideW <- intersect.owin(Xi,w)
  areaXiInside <- area.owin(XiInsideW)
  areaWindow <- area.owin(w)
  
  covProbEstimate <- areaXiInside/areaWindow
  return(covProbEstimate)
}

#' @examples
#' XiOWIN <- heather$coarse
#' windowOWIN <- Frame(heather$coarse)
#' coverageProb <- coveragefrac(XiOWIN,windowOWIN)
#' #
#' #vegetation map of Balcatta park
#' data(balcattapark_coarse)
#' coverageProb <- coveragefrac(balcattapark_coarse$vegmask,balcattapark_coarse$boundary)

#' @keywords spatial nonparametric
#' @seealso \code{\link{exactvarP}} to estimate variance 
