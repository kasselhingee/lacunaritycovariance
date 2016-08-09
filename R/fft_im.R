#' @title Fast Fourier Transform of Images
#' @export fft.im
#' 
#' @description  Uses R's inbuilt FFT function it has addition machinery for returning corrects scales for approximating the continuous fourier transform
#' \deqn{
#'  f(z) = \int g(x) e^{i<x,z>} dx
#' }
#' where \eqn{<x,z>} is the standard spatial dot product: \eqn{<x,z> = \sum_{i=1}^d x_i z_i}.
#' 
#' @param img input image. Can't contain any NA values
#' @param padfactor Optional argument to pad image with zeros and thus increase the spectral resolution (technically interpolation because no extra information is included) of the result.
#' @return an image with correct spectral units and height of the spectral function for approximating a continuous Fourier transform.
#'
fft.im <- function(img, padfactor = c(1,1)) {
  xstep=img$xstep
  ystep=img$ystep
  unitname <- unitname(img)
  M <- as.matrix(img)
  if (any(is.na(M))) {stop("input image, img, contains NA values, cannot fourier transform")}
  stopifnot(is.wholenumber(padfactor[1]) && is.wholenumber(padfactor[2]))
  padfactor=as.integer(padfactor)
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
  
  newunitname <- NULL
  if (!is.null(unitname$singular)){newunitname$singular <- paste0("1/",unitname$singular)}
  if (!is.null(unitname$plural)){newunitname$plural <- paste0("1/",unitname$plural)}
  return(im(fM,
         xcol = xcol,
		 yrow = yrow,
		 unitname= newunitname))
}

#taken from is.integer example page:
is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol