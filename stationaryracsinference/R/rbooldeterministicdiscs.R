#' @title Simulation of Boolean Model of Deterministic Discs
#' @export simulateBooleanDetermDiscs  booldetermdiscs_truecoveragefrac
#' 
#' @description Function for simulating a Boolean model with deterministic discs, and also functions for calculating theoretical properties such as coverage function, covariance, and spectral density.
#' 
#' @param lambda Intensity of the germ process (which is a Poisson point process)
#' @param discr Radius of the discs
#' @param window The window to simulate in (an owin object)
#' @return 
#' \code{simulateBooleanDetermDiscs} returns an owin object containing the locations covered by the Boolean model. The window information is not contained.
#' 
#' \code{booldetermdiscs_truecoveragefrac} returns the true coverage fraction given the intensity and the radius.
#' 
#' 


simulateBooleanDetermDiscs <- function(lambda,discr,window){
  grainlib <- solist(disc(radius=discr))
  bufferdist <- 1.1*discr
  pp <- rpoispp(lambda=lambda,win=dilation(window,bufferdist),nsim=1,drop=TRUE)#lambda from B\"{o}m (2002) - chosen to make coverage probability very close to 0.5
  if (pp$n ==0 ){warning("No points simulated.")
       return(NULL)}
  xibuffer <- placegrainsfromlib(pp,grainlib)
  xi <- intersect.owin(xibuffer,window)
  return(xi)
}

#' @rdname simulateBooleanDetermDiscs
#' 
#' 
booldetermdiscs_truecoveragefrac <- function(lambda, discr){
  return (1-exp(-2*pi*discr^2*lambda))
}

#theoretical set covariance of a disc
# @param r is the radius to calculate set covariance
# @param discr is the radius of disc
setcovdisc <- function(r,discr){
  if (r>=2*discr){setcovariance <- 0}
  else {
    setcovariance <- 2*discr^2*acos(r/(2*discr)) - (r/2)*sqrt(4*discr^2-r^2)
  }
  return(setcovariance)
}

# Isotropic covariance of the Boolean model, given by distance r.
# @param r is the radius to calculate covariance
# @param lambda is the intensity of the germ process (Poisson point process)
# @param discr is the radius of the discs.
isotropiccovarianceDeterministicDiscs <- function(r,lambda,discr){
  expectedsetcovariance <- setcovdisc(r,discr)
  p <- 1-exp(-pi*discr^2*lambda)
  covariance <- 2*p-1+(1-p)^2*exp(lambda*expectedsetcovariance)
  return(covariance)
}
# covariance as a function of vectors given in X, Y columns.
thcovDeterministicDiscs_vec <- function(X,Y,lambda,discr){
  rlist <- sqrt(X^2+Y^2)
  covar <- vector(length(rlist),mode="numeric")
  for (i in 1:length(rlist)){
    covar[i] <- isotropiccovarianceDeterministicDiscs(rlist[i],lambda=lambda,discr=discr)
  }
  return(covar)
}

# get an image of the theoretical covariance
# xrange range of x values
# yrange range of y values
# eps list of distances between samples points in x and y respectively.
thcovDeterministicDiscs <- function(xrange,yrange,eps,lambda,discr){
  xpts <- seq(from = xrange[1], to = xrange[2], by = eps[1])
  ypts <- seq(from = yrange[1], to = yrange[2], by = eps[2])
  mat <- outer(xpts,ypts,FUN="thcovDeterministicDiscs_vec",lambda=lambda,discr=discr) #rows correspond to xstep - just a quirk of outer!
  #reflect out to all corners
  return(im(mat,xcol = xpts,ycol=ypts))
}


#' @section WARNING \code{simulateBooleanDetermDiscs} does not handle the case of an empty realisation very well.
#' This is because the frame of returned value \code{Xi} can be smaller than the simulation window.
#' 
#' 
#' @examples 
#' #Boolean model with discs of radius 10.
#' #The intensity has been chosen such that the true coverage fraction is very close to 0.5.
#' discr <- 10
#' w <- owin(xrange=c(0,100),c(0,100))
#' lambda <- 2.2064E-3
#' xi <- simulateBoolean(lambda,discr,w)
#' plot(xi)
#' plot(w,add=TRUE)