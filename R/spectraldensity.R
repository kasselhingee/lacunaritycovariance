#'Spectral Density of a RACS

#' @description 
#' \code{spectraldensity} estimates the spectral density of a RACS without any kernal smoothing.
#' 
#' @details 
#' Applies FFT to the input image and takes the square of the magnitude to estimate the (non-kernel smoothed) spectral density
#' as described in [B\"{o}hm 2004] (the very first equation).  
#' 
#' 
#' @references Böhm, S., Heinrich, L., Schmidt, V., 2004. Kernel Estimation of the Spectral Density of Stationary Random Closed Sets. Australian & New Zealand Journal of Statistics 46, 41–51. doi:10.1111/j.1467-842X.2004.00310.x
#' 
#' 
#' 
#steal ideas from spatstat's convolve.im

