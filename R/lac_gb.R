#' @title Gliding Box lacunarity from a black and white image
#' @export lacunarityGB
#'
#' @description Calculates the gliding box lacunarity
#'
#' @examples
data(balcattapark_coarse)
img <- balcattapark_coarse$vegmask
bandwidth <- 5 #5 pixel radius?
#create window kernel
  xstep = img$xstep
  ystep = img$ystep
  mat <- matrix(1,ncol=1+2*bandwidth,nrow=1+2*bandwidth)
  kernelfcn <- im(mat,xcol=(-bandwidth:bandwidth)*xstep,yrow=(-bandwidth:bandwidth)*ystep)

  
Algorithm Plan:
+ do FFT (convolve.im) with a uniform kernel.
+ sum over the result in an eroded window
+ change kernel and repeat

+ get lidar data at some point


