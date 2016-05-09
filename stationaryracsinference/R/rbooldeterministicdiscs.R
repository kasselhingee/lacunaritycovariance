#' @title Simulation of Boolean Model of Deterministic Discs
#' @export simulateBooleanDetermDiscs  booldetermdiscs_truecoveragefrac thcovarDeterministicDiscs thspecdensAtOrigin quasithspecdens
#' 
#' @description Function for simulating a Boolean model with deterministic discs, and also functions for calculating theoretical properties such as coverage function, covariance, and spectral density.
#' 
#' @param lambda Intensity of the germ process (which is a Poisson point process)
#' @param discr Radius of the discs
#' @param window The window to simulate in (an owin object)
#' @param seed Optional input (default in NULL). Is an integer passed to \code{\link{base}{set.seed}}. Used to reproduce patterns exactly.

#' @return 
#' \code{simulateBooleanDetermDiscs} returns an owin object containing the locations covered by the Boolean model. The window information is not contained.
#' 
#' \code{booldetermdiscs_truecoveragefrac} returns the true coverage fraction given the intensity and the radius.
#' 
#' \code{thcovarDeterministicDiscs} returns an image of the covariance (which will be isotropic)

#' @section WARNING:
#'  \code{simulateBooleanDetermDiscs} does not handle the case of an empty realisation very well.
#' This is because the frame of returned value \code{Xi} can be smaller than the simulation window.
#' 
#' 
#' @examples 
#' #Boolean model with discs of radius 10.
#' #The intensity has been chosen such that the true coverage fraction is very close to 0.5.
#' discr <- 10
#' w <- owin(xrange=c(0,100),c(0,100))
#' lambda <- 2.2064E-3
#' xi <- simulateBooleanDetermDiscs(lambda,discr,w)
#' plot(xi)
#' plot(w,add=TRUE)
#' 
#' #calculate theoretical values of the model
#' truecoveragefrac <- booldetermdiscs_truecoveragefrac(lambda,discr)
#' truecovariance <- thcovarDeterministicDiscs(
#'                    c(-10,10),c(-10,10),c(0.2,0.2),lambda,discr)
#' thspecdens <- quasithspecdens(lambda,discr)
#' thspecdens_origin <- thspecdensAtOrigin(lambda,discr)
#' thspecdens[round(dim(thspecdens)[2]/2),round(dim(thspecdens)[1]/2)]

simulateBooleanDetermDiscs <- function(lambda,discr,window,seed=NULL){
  grainlib <- solist(disc(radius=discr))
  bufferdist <- 1.1*discr
  
  if (!missing(seed)){set.seed(seed)}
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
  return (1-exp(-pi*discr^2*lambda))
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

#' @rdname simulateBooleanDetermDiscs
#' @param xrange range of x values for \code{thcovarDeterministicDiscs}
#' @param yrange range of y values for \code{thcovarDeterministicDiscs}
#' @param eps list of length 2 of the steps between samples points in x and y respectively for \code{thcovarDeterministicDiscs}.
thcovarDeterministicDiscs <- function(xrange,yrange,eps,lambda,discr){
  xpts <- seq(from = xrange[1], to = xrange[2], by = eps[1])
  ypts <- seq(from = yrange[1], to = yrange[2], by = eps[2])
  mat <- outer(xpts,ypts,FUN="thcovDeterministicDiscs_vec",lambda=lambda,discr=discr) #rows correspond to xstep - just a quirk of outer!
  mat <- t(mat) #now columns correspond to x vals.
  return(im(mat,xcol = xpts, yrow=ypts))
}

#' @describeIn simulateBooleanDetermDiscs  Gives an estimate of the spectral density using the theoretical covariance and FFT
quasithspecdens <- function(lambda,discr){
  xptsLR <- 0:(20*discr)/4
  yptsLR <- 0:(20*discr)/4
  mat <- outer(xptsLR,yptsLR,FUN="thcovDeterministicDiscs_vec",lambda=lambda,discr=discr) #rows correspond to xstep - just a quirk of outer!
  #reflect out to all corners
  mat <- mat[,c((ncol(mat)):2,1:ncol(mat))]
  mat <- mat[c((nrow(mat)):2,1:nrow(mat)),]
  #theoretical p value 
  p <- 1-exp(-pi*discr^2*lambda)
  
  M <- mat-p^2
  nr <- nrow(M)
  nc <- ncol(M)
  #pad with lots of 0!
  thcovpad <- matrix(0, ncol=8*nc, nrow=8*nr)
  thcovpad[1:nr, 1:nc] <- M
  scalefactorX <- (xptsLR[2]-xptsLR[1])
  scalefactorY <- (yptsLR[2]-yptsLR[1]) 
  specdens <- scalefactorX*scalefactorY*fft(thcovpad)
  specdens <- abs(specdens)
  nr <- nrow(specdens) #can use these because specdens has same dimensions as start
  nc <- ncol(specdens)
  if (nr %% 2 == 0){
    specdens <- specdens[ ((-nr/2):(nr/2)) %% (nr) + 1,]
    yrow <- ((-nr/2):(nr/2)) * 2*pi/(scalefactorY*nr)
  } else {
    specdens <- specdens[ ((-(nr-1)/2):((nr-1)/2)) %% (nr) + 1,]
    yrow <- ((-(nr-1)/2):((nr-1)/2)) * 2*pi/(scalefactorY*nr)
  } 
  if (nc %% 2 == 0){
    specdens <- specdens[, ((-nc/2):(nc/2)) %% (nc) + 1]
    xcol <- ((-nc/2):(nc/2)) * 2*pi/(scalefactorX*nc)
  } else {
    specdens <- specdens[, ((-(nc-1)/2):(nc/2)) %% (nc) + 1]  
    xcol <- ((-(nc-1)/2):(nc/2)) * 2*pi/(scalefactorX*nc)
  }
  
  specdens <- im(specdens,xcol = xcol, yrow = yrow)
  return(specdens)
}

#' @describeIn simulateBooleanDetermDiscs Calculates spectral density at \eqn{o} using the theoretical covariance.
#'  It is the integral of the covariance - (coverage fraction)^2 over all space.
#now calculate spectral density at origin using integral of covariance -p^2
thspecdensAtOrigin <- function(lambda,discr){
  xptsLR <- 0:(80*discr)/20
  yptsLR <- 0:(80*discr)/20
  mat <- outer(xptsLR,yptsLR,FUN="thcovDeterministicDiscs_vec",lambda=lambda,discr=discr) #rows correspond to xstep - just a quirk of outer!
  #reflect out to all corners
  mat <- mat[,c((ncol(mat)):2,1:ncol(mat))]
  mat <- mat[c((nrow(mat)):2,1:nrow(mat)),]
  #theoretical p value 
  p <- 1-exp(-pi*discr^2*lambda)
  
  M <- mat-p^2
  nr <- nrow(M)
  nc <- ncol(M)
  scalefactorX <- (xptsLR[2]-xptsLR[1])
  scalefactorY <- (yptsLR[2]-yptsLR[1]) 
  return(sum(M)*scalefactorY*scalefactorX)
}



 