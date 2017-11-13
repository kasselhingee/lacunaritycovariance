#' @title Variance Estimates for Observed Area
#' @export varofobsarea.v1
#' @description Estimates the variance of the area of a cover type observed in a thematic map created using a fallible classifier from remote sensing.
#' @author{Kassel Hingee}



#' @param xi An observation of the RACS of interest.
#' @param obswin 

#' @examples
#' xi <- as.im(heather$coarse, na.replace = 0)
#' obswin <- Frame(xi)
#' corrrad <- 10 #10 pixels
#' corrstepheight <- 0.9
#' ommrate <- 0.05
#' commrate <- 0.01
#' xi.sum <- sum(xi)
#' exparea <- xi.sum - xi.sum * ommrate + (sum(1 - xi) * commrate)
#' varguess <- varofobsarea.v1(xi, obswin, corrrad, corrstepheight, ommrate, commrate)


#' @details 
#' To install OpenImageR had to install libtiff5-dev on my ubuntu machine
#' @references 

varofobsarea.v1 <- function(xi, obswin, corrrad, corrstepheight, ommrate, commrate){
  xi[complement.owin(obswin, frame = Frame(xi))] <- 0
  #radius filter of the cover type of interest
  xiconvsum <- convandintersectsum(xi, corrrad)
  varfromomm <- ommrate * (1 - ommrate) * ( (1 - corrstepheight) * sum(xi) + corrstepheight * xiconvsum)
  
  #radius filter of the alternate cover type
  xic <- 1-xi
  xic[complement.owin(obswin, frame = Frame(xi))] <- 0  
  xicconvsum <- convandintersectsum(xic, corrrad)
  varfromcomm <- commrate * (1 - commrate) * ( (1 - corrstepheight) * sum(xic) + corrstepheight * xicconvsum)
  
  return(varfromomm + varfromcomm)
}

convandintersectsum <- function(xi, corrrad){
  # create disc kernel
  kernel.owin <- disc(radius = corrrad, mask = TRUE, eps = 1)#make everything in pixel units
  kernel.m <- as.matrix(kernel.owin)

  xi.m <- as.matrix(xi)
  
  
  # apply radius filter to xi
  convolved <- OpenImageR::convolution(xi.m, kernel.m, mode = "same")
  
  # sum over filter results by points in xi
  convolved[!xi.m] <- 0
  return(sum(convolved))
}
