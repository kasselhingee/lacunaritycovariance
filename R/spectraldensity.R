#'Spectral Density of a RACS

#' @description 
#' \code{spectraldensity} estimates the spectral density of a RACS without any kernal smoothing.
#' 
#' @details 
#' Applies FFT to the input image and takes the square of the magnitude to estimate the (non-kernel smoothed) spectral density
#' as described in [B\"{o}hm 2004] (the very first equation or equivalently \eqn{\hat{h}(x)} of equation 4.6 in the 2002 report).  
#' 
#' 
#' @references Böhm, S., Heinrich, L., Schmidt, V., 2004. Kernel Estimation of the Spectral Density of Stationary Random Closed Sets. Australian & New Zealand Journal of Statistics 46, 41–51. doi:10.1111/j.1467-842X.2004.00310.x
#' 
#' 
#' @param X A rectangular observation window. NA's are assumed to mean outside \eqn{\Xi} rather than missing data or anything (**I haven't assesed the theory/computations for non-rectangular windows but its probably the same)
#steal ideas from spatstat's convolve.im
spectraldensity <- function(X){
  stopifnot(is.im(X))
  M <- X$v
  M[is.na(M)] <- 0 #since the function that we wish to transform is an indicator of both inside window, and inside xi. Its ok to set all NAs to 0
  M <- M-p
  fM <- fft(M)
  nr <- nrow(M)
  nc <- ncol(M)
  areaM <- nr * nc #because theory uses rectangular windows, I'm going to assume a rectangular window to - maybe improve on this later
  specdens <- (Re(fM)^2+Im(fM)^2)/(areaM^3) #divide by areaM^2 to put integrals on correct scale. Divide by areaM again to get to spectral density
  #currently specdens[i,j] corresponds to a spectral location of 
  #     y = ((i-1) mod numrow)/ystep, x = ((j-1) mod numcol)/xstep
  # Rearrange this periodic function so that 
  # the origin of translations (0,0) is at matrix position (nr/2,nc/2) or close depending on whether nr is even or not
  # NB this could introduce an extra row and column
  if (nr %% 2 == 0){
    specdens <- specdens[ ((-nr/2):(nr/2)) %% (nr) + 1,]
    yrow <- ((-nr/2):(nr/2)) * 1/X$ystep
  } else {
    specdens <- specdens[ ((-(nr-1)/2):((nr-1)/2)) %% (nr) + 1,]
    yrow <- ((-(nr-1)/2):((nr-1)/2)) * 1/X$ystep
  } 
  if (nc %% 2 == 0){
    specdens <- specdens[, ((-nc/2):(nc/2)) %% (nc) + 1]
    xcol <- ((-nc/2):(nc/2)) * 1/X$xstep
  } else {
    specdens <- specdens[, ((-(nc-1)/2):(nc/2)) %% (nc) + 1]  
    xcol <- ((-(nc-1)/2):(nc/2)) * 1/X$xstep
  }

  specdens <- im(specdens,xcol = xcol, yrow = yrow)
  
  return(specdens)
  
}
