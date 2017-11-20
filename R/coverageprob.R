#' @title Estimate the coverage probability of a stationary RACS
#' @author Kassel Hingee 
#' @import spatstat
#' @export coveragefrac coverageprob
#' 
#' @description 
#' The coverage probability of a stationary RACS is also known as the coverage fraction and is the probability that an arbitrary point is covered by the RACS.
#' Given a realisation, \code{xi}, in window \eqn{W} of a stationary RACS this function computes the fraction of \eqn{W} covered by \code{xi}, which is an estimate of the coverage probability.
#' This estimate relies on the window \eqn{W} being large compared to the spatial dependence of the RACS.
#' See [1, section 6.4.2] for more theoretical details.
#' 
#' @references [1] Chiu, S.N., Stoyan, D., Kendall, W.S. and Mecke, J. (2013) Stochastic Geometry and Its Applications, 3rd ed. Chichester, United Kingdom: John Wiley & Sons.
 
#' 
#' @param xi An observation of a RACS in \pkg{spatstat}'s \code{owin} or \code{im} format.
#' @param obswin The window of observation (not necessarily rectangular) also in \code{owin} format.
#' @return An estimate of the coverage probability
#' @details
#' If \code{xi} is in \code{im} format then \code{xi} must be an image of 1s, 0s and NAs
#'  representing inside the set, outside the set and outside the observation window respectively.
#'  \code{coverageprob} will not accept a \code{obswin} argument if \code{xi} is in \code{im} format.
#' 
#' @examples
#' xi <- heather$coarse
#' obswindow <- Frame(heather$coarse)
#' cp <- coverageprob(xi, obswindow)

#' @keywords spatial nonparametric
coverageprob <- function(xi, obswin = NULL){
  if (is.im(xi)){
    stopifnot(is.null(obswin))
    coverprobest <- sum(xi) / sum(is.finite(xi$v))
    return(coverprobest)
  }
  stopifnot(is.owin(xi))
  stopifnot(is.owin(obswin))
  xinw <- intersect.owin(xi, obswin)
  xiinw_area <- area.owin(xinw)
  w_area <- area.owin(obswin)

  coverprobest <- xiinw_area / w_area
  return(coverprobest)
}

#' @rdname coverageprob 
coveragefrac <- coverageprob

