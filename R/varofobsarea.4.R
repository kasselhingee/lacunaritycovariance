#' @title Variance Estimates for Observed Area
#' @export sae.v4.mean sae.v4.var
#' @description Estimates the variance of the area of a cover type observed in a thematic map created using a fallible classifier from remote sensing.
#' @author{Kassel Hingee}



#' @param xi An observation of the RACS of interest in owin form.
#' @param obswin Observation window in owin format
#' @param corrrad Radius of the step function in the correlation (in the same units as xi)
#' @param corrstepheight Height of the step in the correlation
#' @param p21 Probability that a location fallibly classified into the cover of interest is really outside the cover of interest
#' @param p12 Probability that a location fallibly classified into the non-interesting class is really the class of interest.


#' @examples
#' xi <- heather$coarse
#' obswin <- Frame(xi)
#' erosionrad <- 0.1 #in units of image
#' exparea <- sae.v4.mean(xi, obswin)
#' varguess <- sae.v4.var(xi, obswin, erosionrad)
#' stdofareaest <- sqrt(varguess)


#' @details 
#' To install OpenImageR had to install libtiff5-dev on my ubuntu machine
 
sae.v4.mean <- function(xi, obswin){
  xi <- intersect.owin(xi, obswin)
  xi.area <- area.owin(xi)
  window.area <- area.owin(obswin)
  return(xi.area)
}

#' @describeIn sae.v4.mean A `variance' estimated by adhoc selection of erosion and dilation distances corresponding to 2*standard deviation
sae.v4.var <- function(xi, obswin, erosionrad){
  xi <- intersect.owin(xi, obswin)
  xi.c <- setminus.owin(obswin, xi)
  xi.erode <- setminus.owin(xi, dilation(xi.c, erosionrad))

  xi.dilate <- dilation(xi, erosionrad)
  xi.dilate <- intersect.owin(xi.dilate, obswin)

  xi.area <- area.owin(xi)
  twosigma <- max(xi.area - area.owin(xi.erode), area.owin(xi.dilate) - xi.area)
  
  return((twosigma/2)^2)
}
