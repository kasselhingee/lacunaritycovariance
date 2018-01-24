#' @title Variance Estimates for Observed Area
#' @export varofobsarea.v3 expectedarea
#' @description Estimates the variance of the area of a cover type observed in a thematic map created using a fallible classifier from remote sensing.
#' @author{Kassel Hingee}



#' @param xi An observation of the RACS of interest in 1, 0, or NA valued pixels.
#' Pixels must be square. Must be a spatstat im object.
#' @param obswin Observation window
#' @param corrrad Radius of the step function in the correlation (in the same units as xi)
#' @param corrstepheight Height of the step in the correlation
#' @param p01 Probability of an randomly chosen fallibly-classified
#'  tree pixel is not in the true tree canopy
#' @param p10 Probability that a falliblty classified non-tree pixel is 
#' in the tree canopy.

#' @examples
#' xi <- as.im(heather$coarse, na.replace = 0)
#' obswin <- Frame(xi)
#' corrrad <- 5 #in units of image
#' corrstepheight <- 0.9
#' p01 <- 0.05
#' p10 <- 0.01
#' xi.sum <- sum(xi)
#' exparea <- expectedarea(xi, obswin, p01 = p01, p10 = p10)
#' varguess <- varofobsarea.v3(xi, obswin, corrrad, corrstepheight, p01, p10)
#' stdofareaest <- sqrt(varguess)


#' @details 
#' To install OpenImageR had to install libtiff5-dev on my ubuntu machine
 
expectedarea <- function(xi, obswin, p01=NA, p10=NA){
  xi <- xi[as.rectangle(obswin), drop=TRUE]
  xi[complement.owin(obswin)] <- 0
  xic <- 1-xi
  xic[complement.owin(obswin)] <- 0  
  
  return((sum(xi) * (1 - p01) + sum(xic) * p10) * xi$xstep * xi$ystep)
}

#' @describeIn expectedarea The variance of the expected area using Small Area Estimation Method - Version 3. 
#' This methods assumes that the omission and comission errors are independent and different processes and that the correlation between errors (within each of theses processes) is a step function with radius \code{corrrad}.
varofobsarea.v3 <- function(xi, obswin, corrrad, corrstepheight, p01, p10){
  xi <- xi[as.rectangle(obswin), drop=TRUE]
  xi[complement.owin(obswin)] <- 0
  #radius filter of the cover type of interest
  xiconvsum <- convandintersectsum(xi, corrrad)
  varfromomm <- p01 * (1 - p01) * ( (1 - corrstepheight) * sum(xi) + corrstepheight * xiconvsum)
  
  #radius filter of the alternate cover type
  xic <- 1-xi
  xic[complement.owin(obswin)] <- 0  
  xicconvsum <- convandintersectsum(xic, corrrad)
  varfromcomm <- p10 * (1 - p10) * ( (1 - corrstepheight) * sum(xic) + corrstepheight * xicconvsum)
  
  return((varfromomm + varfromcomm) * xi$xstep^2 * xi$ystep^2)
}

convandintersectsum <- function(xi, corrrad){
  # create disc kernel
  kernel.owin <- disc(radius = corrrad, mask = TRUE, eps = c(xi$xstep, xi$ystep)) #make disc at same resolution as xi
  kernel.m <- as.matrix(kernel.owin)

  xi.m <- as.matrix(xi)
  
  
  # apply radius filter to xi
  convolved <- OpenImageR::convolution(xi.m, kernel.m, mode = "same")
  
  # sum over filter results by points in xi
  convolved[!xi.m] <- 0
  return(sum(convolved))
}
