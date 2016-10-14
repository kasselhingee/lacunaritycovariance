#' @title Gliding Box lacunarity from a black and white image
#' @export lacunarityGB
#'
#' @description Calculates the gliding box lacunarity
#'
#' @examples
library(stationaryracsinference)
data(balcattapark_coarse)
img <- balcattapark_coarse$vegmask
img <- as.im(balcattapark_coarse$vegmask)
bandwidth <- 5 #5 pixel radius?
#create window kernel
  xstep = img$xstep
  ystep = img$ystep
  mat <- matrix(1/((1+2*bandwidth)^2),ncol=1+2*bandwidth,nrow=1+2*bandwidth)
  kernelfcn <- im(mat,xcol=(-bandwidth:bandwidth)*xstep,yrow=(-bandwidth:bandwidth)*ystep)

  convolvedimg <- convolve.im(img,kernelfcn)
  dim(convolvedimg)
  dim(img)

#lacunarity if you let the box centeres be everywhere
  smA <- mean(convolvedimg) #sample mean
  ss2A <- mean(convolvedimg^2) #biased sample second moment
  lacA <- ss2A/(smA^2) -1
lacA
  convolvedimgRS <-  as.im(convolvedimg,W=erosion(Window(img),bandwidth))
  Window(convolvedimgRS)
  Window(img)
  smB <- mean(convolvedimgRS) #sample mean
  ss2B <- mean(convolvedimgRS^2) #biased sample second moment
  lacB <- ss2B/(smB^2) -1
lacB
plot(img)
plot(convolvedimg)
plot(convolvedimgRS)
  
Algorithm Plan:
+ do FFT (convolve.im) with a uniform kernel.
+ sum over the result in an eroded window
+ change kernel and repeat

+ get lidar data at some point


