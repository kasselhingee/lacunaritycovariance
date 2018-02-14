#' @title Variance Estimates for Observed Area
#' @export sae.v3.mean sae.v3.var
#' @description Estimates the variance of the area of a cover type observed in a thematic map created using a fallible classifier from remote sensing.
#' @author{Kassel Hingee}



#' @param xi An observation of the RACS of interest in 1, 0, or NA valued pixels.
#' Pixels must be square. Must be a spatstat im object.
#' @param obswin Observation window
#' @param corrrad Radius of the step function in the correlation (in the same units as xi)
#' @param corrstepheight Height of the step in the correlation
#' @param p21 Probability that a location fallibly classified into the cover of interest is really outside the cover of interest
#' @param p12 Probability that a location fallibly classified into the non-interesting class is really the class of interest.


#' @examples
#' xi <- as.im(heather$coarse, na.replace = 0)
#' obswin <- Frame(xi)
#' corrrad <- 5 #in units of image
#' corrstepheight <- 0.9
#' p21 <- 0.05
#' p12 <- 0.01
#' xi.sum <- sum(xi)
#' exparea <- sae.v3.mean(xi, obswin, p21 = p21, p12 = p12)
#' varguess <- sae.v3.var(xi, obswin, corrrad, corrstepheight, p21, p12)
#' stdofareaest <- sqrt(varguess)


#' @details 
#' To install OpenImageR had to install libtiff5-dev on my ubuntu machine
 
sae.v3.mean <- function(xi, obswin, p21=NA, p12=NA){
  xi <- xi[as.rectangle(obswin), drop=TRUE]
  xi[complement.owin(obswin, frame = Frame(dilation.owin(obswin,1)))] <- 0
  xic <- 1-xi
  xic[complement.owin(obswin, frame = Frame(dilation.owin(obswin,1)))] <- 0  
  
  return((sum(xi) * (1 - p21) + sum(xic) * p12) * xi$xstep * xi$ystep)
}

#' @describeIn sae.v3.mean The variance of the expected area using Small Area Estimation Method - Version 3. 
#' This methods assumes that the omission and comission errors are independent and different processes and that the correlation between errors (within each of theses processes) is a step function with radius \code{corrrad}.
sae.v3.var <- function(xi, obswin, corrrad, corrstepheight, p21, p12){
  xi <- xi[as.rectangle(obswin), drop=TRUE]
  xi[complement.owin(obswin, frame = Frame(dilation.owin(obswin,1)))] <- 0
  #radius filter of the cover type of interest
  xiconvsum <- convandintersectsum(xi, corrrad)
  varfromcomm <- p21 * (1 - p21) * ( (1 - corrstepheight) * sum(xi) + corrstepheight * xiconvsum)
  
  #radius filter of the alternate cover type
  xic <- 1-xi
  xic[complement.owin(obswin, frame = Frame(dilation.owin(obswin,1)))] <- 0  
  xicconvsum <- convandintersectsum(xic, corrrad)
  varfromomm <- p12 * (1 - p12) * ( (1 - corrstepheight) * sum(xic) + corrstepheight * xicconvsum)
  
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


