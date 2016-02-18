#' @title Kernel Smoothing Tools
#' @export kernelsmooth
#' 
#' @description Function for kernel smoothing spatstat images. This is mostly for smoothing spectral density estimates.
#' Only kernel currently supported is Epanechnikov.
#' 
#' @details Uses FFT to do the calculation.
#' @param bandwidth Bandwidth for the chosen kernel (in reciprocal of the units of \code{im})
#' @param kernel Specifies the kernel to use. Currently only Epanechnikov is supported (and is the default).
#' @param im an image for smoothing.
#' @return A smoothed image in spatstat \code{im} format
kernelsmooth <- function(im,bandwidth,kernel="Epanechnikov"){
  stopifnot(is.im(im))
  xstep = im$xstep
  ystep = im$ystep
  if (kernel == "Epanechnikov") {supportwidth = bandwidth}
  else {
    warning("In spectral density smoothing support width of kernel is unknown, defaulting to 3x bandwidth")
    supportwidth = 3*bandwidth
  }
  X <- seq(0,supportwidth*1.5+xstep,by=xstep) #much larger than support width to avoid boundary issues?
  Y <- seq(0,supportwidth*1.5+ystep,by=ystep)
  if (kernel == "Epanechnikov"){
    mat <- outer(X/bandwidth,Y/bandwidth,FUN="EpanechnikovFcn") #rows correspond to xstep - just a quirk of outer!
    mat <- mat/(bandwidth^2) #to account for the scaling of the kernel - so that it all adds to 1
    mat <- t(mat) #columns correspond to changes in X, rows correspond to changes in Y!
    #reflect out to all corners
    mat <- mat[,c((ncol(mat)):2,1:ncol(mat))]
    mat <- mat[c((nrow(mat)):2,1:nrow(mat)),]
    kernelfcn <- im(mat,xcol=c(-X[length(X):2],X),yrow=c(-Y[length(Y):2],Y))
  }
  else {
    stop("Given kernel is not supported for spectral density smoothing")
  }
  #apply kernel using convolve.im
  smim <- convolve.im(im,kernelfcn)
  smim <- smim[Frame(im)]
  return(smim)
}

#EpanechnikovFcn 
EpanechnikovFcn <- function(X,Y){#WARNING: operates in 2D only on a vector of things 
  stopifnot(length(X)==length(Y))
  result <- vector(length=length(X),mode="numeric")
  sz <- sqrt((X*X)+(Y*Y))
  result[sz>1] <- 0
  result[sz<=1] <- (2/pi)*(1-sz[sz<=1]^2)
  return(result)
}

#' @examples 
#' #smooth a step image
#' stepim <- im(matrix(0,nrow=100,ncol=100),xrange=c(-1,1),yrange=c(-1,1))
#' stepim[owin(xrange=c(-0.5,0.5),yrange=c(-0.5,0.5))] <- 1
#' plot(stepim,axes=TRUE)
#' 
#' smstepim <- kernelsmooth(stepim,0.5)
#' plot(smstepim,axes=TRUE)