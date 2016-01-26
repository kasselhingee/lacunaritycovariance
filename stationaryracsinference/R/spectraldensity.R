#'Spectral Density of a RACS
#' @export unsmoothedspectraldensity spectraldensity
#' 
#' @description 
#' \code{unsmoothedspectraldensity} estimates the spectral density of a RACS without any kernal smoothing.
#' \code{spectraldensity} smoothes the estimate and is my recommended function due to known asymptotoic properties [1].
#' 
#' @details 
#' Applies FFT to the input image and takes the square of the magnitude to estimate the (non-kernel smoothed) spectral density
#' as described in [B\"{o}hm 2004] (the very first equation or equivalently \eqn{\hat{h}(x)} of equation 4.6 in the 2002 report).  
#' 
#' \eqn{\hat{h}(z) := \frac{1}{|w|}\left|\int_{R^2}1_w(x) (1_\Xi(x)-p) dx\right|}
#' 
#' This is what \code{unsmoothedspectraldensity} returns.
#'  \code{spectraldensity} then kernel smooths \eqn{\hat{h}} before returning a spectral density estimate.
#' 
#' @references B\"{o}hm, S., Heinrich, L., Schmidt, V., 2004. Kernel Estimation of the Spectral Density of Stationary Random Closed Sets. Australian & New Zealand Journal of Statistics 46, 41--51. doi:10.1111/j.1467-842X.2004.00310.x
#' 
#' 
#' @param Xi A rectangular observation of \eqn{Xi}. NA's are assumed to mean outside \eqn{\Xi} rather than missing data. Xi must be pixel mask owin object. (**I haven't assesed the theory/computations for non-rectangular windows but its probably the same)
#' @param w Observation window (must be rectangular for now)
#' @param bandwidth Bandwidth for the chosen kernel (in reciprocal of the units of \code{Xi})
#' @param kernel Specifies the kernel to use. Currently only Epanechnikov is supported (and is the default).
#' @param ... Arguments passed to \code{as.mask} to convert observation into a pixel image. This is useful for specifying the resolution of the raster derived from \code{Xi} (otherwise the raster will always have 128x128 pixels).
#' 
#' 
#' 
#' @section WARNING: WARNING RESULTS OF THIS FUNCTION HAVENT BEEN TESTED - THEY COULD BE WILDLY WRONG
#' 
#' Also final result should have odd pixel dimensions (so that the spectral density at origin is easy to pull out), but convolve.im returns even dimensions, I don't know why yet.
#' 

#' @examples 
#' #apply to heather
#' xi <- heather$coarse
#' specdens <- spectraldensity(xi,Frame(xi),20)
#' plot(specdens)
#' 
#' #See what unsmoothed version looks like
#' unsmspecdens <- unsmoothedspectraldensity(xi,Frame(xi))
#' plot(unsmspecdens)
#' plot(solist(smoothed = specdens,unsmoothed = unsmspecdens),
#'        main = "Spectral Density Estimates",
#'        equal.scales=TRUE)
#' 



spectraldensity <- function(Xi,w,bandwidth,kernel="Epanechnikov",...){
  specdens <- unsmoothedspectraldensity(Xi,w,...)
  xstep = specdens$xstep
  ystep = specdens$ystep
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
  smspecdens <- convolve.im(specdens,kernelfcn)
  smspecdens <- smspecdens[Frame(specdens)]
  return(smspecdens)
}


#steal ideas from spatstat's convolve.im
#' @rdname spectraldensity
#' @param suffspecres Optional. A desired spectral resolution, the output resolution will be smaller than \code{suffspecres}
#'  If it is not provided then resolution will be 2*pi/(size of spatial extent of window), for example a window 500m x 500m will result in spectral resolution of 2*pi/500 in both X and Y directions.
#'  Can be of length 1 (applies to all dimensions), or length 2 (a differen desired resolution for each dimension).
unsmoothedspectraldensity <- function(Xi,w,suffspecres=NULL,...){
  stopifnot(is.owin(Xi))
  stopifnot(is.mask(Xi))
  stopifnot(is.rectangle(w)) #because theory uses rectangular windows, I'm going to assume a rectangular window to - maybe improve on this later
  Xi <- intersect.owin(Xi,w)
  p <- covpest(Xi,w)
  xstep=Xi$xstep
  ystep=Xi$ystep
  M <- as.matrix(Xi)
  M[is.na(M)] <- 0 #since the function that we wish to transform is an indicator of inside window AND inside xi. Its ok to set all NAs to 0
  M <- M-p
  #pad if necessary #DFT approximation makes the assumption that function is 0 outside window anyway
  if ((!is.null(suffspecres))
    &&
    ( any(suffspecres < 2*pi/(c(xstep*ncol(M),ystep*nrow(M)))) ))
    {
      #calculate padding required in each direction: kX such that suffspecresX > 2*pi/(xstep*nrow(M)*kX) (so output BETTER than suffspecres)
      padfactor <- ceiling(2*pi/(c(xstep*ncol(M),ystep*nrow(M))*suffspecres))
      Mpad <- matrix(0,nrow=nrow(M)*padfactor[2],ncol=ncol(M)*padfactor[1])
      Mpad[1:nrow(M),1:ncol(M)] <- M
      fM <- xstep*ystep*fft(Mpad) #xstep*ystep = scale to approximate Fourier transform
  }
  else {
      fM <- xstep*ystep*fft(M) #xstep*ystep = scale to approximate Fourier transform
  }
  areaM <- xstep*ncol(M) * ystep*nrow(M)
  specdens <- (Re(fM)^2+Im(fM)^2)/(areaM) #divide by areaM to get spectral density (formula in Bohm)
  #currently specdens[i,j] corresponds to a spectral location of 
  #     y = 2pi*((i-1) mod numrow)/(length*ystep), x = 2pi*((j-1) mod numcol)/(length*xstep)
  # Rearrange this periodic function so that 
  # the origin of translations (0,0) is at matrix position (nr/2,nc/2) or close depending on whether nr is even or not
  # NB this could introduce an extra row and column
  nr <- nrow(fM)
  nc <- ncol(fM)
  if (nr %% 2 == 0){
    specdens <- specdens[ ((-nr/2):(nr/2)) %% (nr) + 1,]
    yrow <- ((-nr/2):(nr/2)) * 2*pi/(ystep*nr)
  } else {
    specdens <- specdens[ ((-(nr-1)/2):((nr-1)/2)) %% (nr) + 1,]
    yrow <- ((-(nr-1)/2):((nr-1)/2)) * 2*pi/(ystep*nr)
  } 
  if (nc %% 2 == 0){
    specdens <- specdens[, ((-nc/2):(nc/2)) %% (nc) + 1]
    xcol <- ((-nc/2):(nc/2)) * 2*pi/(xstep*nr)
  } else {
    specdens <- specdens[, ((-(nc-1)/2):(nc/2)) %% (nc) + 1]  
    xcol <- ((-(nc-1)/2):(nc/2)) * 2*pi/(xstep*nr)
  }

  specdens <- im(specdens,xcol = xcol, yrow = yrow)
  
  return(specdens)
  
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
