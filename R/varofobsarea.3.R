#' @title Variance Estimates for Observed Area
#' @export sae.v3.mean sae.v3.var sae.v3.wsu.mean sae.v3.wsu.var
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
  xi.area <- area.owin(xi)
  window.area <- area.owin(obswin)
  return(xi.area - ( xi.area * p21 ) + ( (window.area - xi.area) * p12))
}

#' @describeIn  sae.v3.mean Expected area including sampling uncertainty
sae.v3.wsu.mean <- function(xi, obswin, n11, n21, n12, n22){
  xi.area <- area.owin(xi)
  window.area <- area.owin(obswin)
  return(xi.area * ahat11(n11, n21) + ( (window.area - xi.area) * ahat12(n12, n22)))
}

#' @describeIn sae.v3.mean The variance of the expected area using Small Area Estimation Method - Version 3. 
#' This methods assumes that the omission and comission errors are independent and different processes and that the correlation between errors (within each of theses processes) is a step function with radius \code{corrrad}.
sae.v3.var <- function(xi, obswin, corrrad, corrstepheight, p21, p12){
  xi.im <- as.im(xi, na.replace = 0) #warning is 0 outside obswin too!
  xi.im <- xi.im[as.rectangle(obswin), drop=TRUE]
  #radius filter of the cover type of interest
  xiconvsum <- convandintersectsum(xi.im, corrrad)
  varfromcomm <- p21 * (1 - p21) * ( (1 - corrstepheight) * sum(xi.im) + corrstepheight * xiconvsum) * xi.im$xstep^2 * xi.im$ystep^2
  
  #radius filter of the alternate cover type
  xi.c <- setminus.owin(obswin, xi) #I should do this instead of complement.owin.inwindow
  xi.c.im <- as.im(xi.c, na.replace = 0)
  xi.c.im <- xi.c.im[as.rectangle(obswin), drop=TRUE]
  xicconvsum <- convandintersectsum(xi.c.im, corrrad)
  varfromomm <- p12 * (1 - p12) * ( (1 - corrstepheight) * sum(xi.c.im) + corrstepheight * xicconvsum) * xi.c.im$xstep^2 * xi.c.im$ystep^2
  
  return(varfromomm + varfromcomm)
}

#' @describeIn sae.v3.mean The variance of the expected area using Small Area Estimation Method - Version 3 and sampling uncertainty in conditional confusion matrix 
#' This methods assumes that the omission and comission errors are independent and different processes and that the correlation between errors (within each of theses processes) is a step function with radius \code{corrrad}.
sae.v3.wsu.var <- function(xi, obswin, corrrad, corrstepheight, n11, n21, n12, n22){
  xi.im <- as.im(xi, na.replace = 0) #warning is 0 outside obswin too!
  warning("owin converted to image using default pixel amounts")
  xi.im <- xi.im[as.rectangle(obswin), drop=TRUE]
  #radius filter of the cover type of interest
  xiconvsum <- convandintersectsum(xi.im, corrrad)
  a11 <- ahat11(n11, n21)
  a12 <- ahat12(n12, n22)
  expectvarfromcomm <- (1 - a11) * a11 * ( (1 - corrstepheight) * sum(xi.im) + corrstepheight * xiconvsum) * xi.im$xstep^2 * xi.im$ystep^2

  #radius filter of the alternate cover type
  xi.c <- setminus.owin(obswin, xi) #I should do this instead of complement.owin.inwindow
  xi.c.im <- as.im(xi.c, na.replace = 0)
  warning("owin converted to image using default pixel amounts")
  xi.c.im <- xi.c.im[as.rectangle(obswin), drop=TRUE]
  xicconvsum <- convandintersectsum(xi.c.im, corrrad)
  expectvarfromomm <- a12 * (1 - a12) * ( (1 - corrstepheight) * sum(xi.c.im) + corrstepheight * xicconvsum) * xi.c.im$xstep^2 * xi.c.im$ystep^2
  
  xi.area <- area.owin(xi)
  window.area <- area.owin(obswin)
  varofexpect <- xi.area^2 * ssq11(n11, n21) + (window.area - xi.area)^2 * ssq12(n12, n22)
  
  return(expectvarfromomm + expectvarfromcomm + varofexpect)
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


