#' @title Fast Fourier Transform of Images
#' 
#' @description  Uses R's inbuilt FFT function it has addition machinery for returning corrects scales for approximating the continuous fourier transform
#' \deqn{
#'  f(z) = \int g(x) e^{i<x,z>} \dx
#' }
#' where \eqn{<x,z>} is the standard spatial dot product: \eqn{<x,z> = \sum_{i=1}^d x_i z_i}.
#' 
#' @param image input image. Can't contain any NA values
#' @param padfactor Optional argument to pad image with zeros and thus increase the spectral resolution (technically interpolation because no extra information is included) of the result.
#' @return an image with correct spectral units.
fft.im <- function(image, padfactor = c(1,1)) {
  xstep=Xi$xstep
  ystep=Xi$ystep
  M <- as.matrix(Xi)
  if (any(is.na(M))) {stop("image contains NA values, cannot fourier transform")}
  stopifnot(is.integer(padfactor[1]) && is.integer(padfactor[2]))
  if (any(padfactor > 1)){
    Mpad <- matrix(0,nrow=nrow(M)*padfactor[2],ncol=ncol(M)*padfactor[1])
    Mpad[1:nrow(M),1:ncol(M)] <- M
    fM <- xstep*ystep*fft(Mpad) #xstep*ystep = scale to approximate Fourier transform
  }
  else {
    fM <- xstep*ystep*fft(M) #xstep*ystep = scale to approximate Fourier transform
  }
  
  #convert to an image with correct dimensions, and centred on spectral value of 0
  nr <- nrow(fM) #can use these because fM has same dimensions as input matrix start
  nc <- ncol(fM)
  if (nr %% 2 == 0){
    fM <- fM[ ((-nr/2):(nr/2)) %% (nr) + 1,]
    yrow <- ((-nr/2):(nr/2)) * 2*pi/(ystep*nr)
  } else {
    fM <- fM[ ((-(nr-1)/2):((nr-1)/2)) %% (nr) + 1,]
    yrow <- ((-(nr-1)/2):((nr-1)/2)) * 2*pi/(ystep*nr)
  } 
  if (nc %% 2 == 0){
    fM <- fM[, ((-nc/2):(nc/2)) %% (nc) + 1]
    xcol <- ((-nc/2):(nc/2)) * 2*pi/(xstep*nc)
  } else {
    fM <- fM[, ((-(nc-1)/2):(nc/2)) %% (nc) + 1]  
    xcol <- ((-(nc-1)/2):(nc/2)) * 2*pi/(xstep*nc)
  }
  
  return(im(fM,xcol = xcol, yrow = yrow))
}