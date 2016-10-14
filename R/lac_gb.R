#' @title Gliding Box lacunarity from a black and white image
#' @export lacunarityGB
#'
#' @description Calculates the gliding box lacunarity
#' @details Calculates the gliding box lacunarity for a given range of box sizes (`radius').
#' The centres are given the pixels of \code{img}.
#' @references **Algorithm described in Allain1991ch**
#'
#' @return An fv object with items for no edge correction and reduced sample border correction
#' @param img An image of 0's and 1's.
#' @param b Bandwidth in the dimensions of \code{img}
#' @examples
#' 
library(stationaryracsinference)
data(balcattapark_coarse)
img <- balcattapark_coarse$vegmask
img <- as.im(balcattapark_coarse$vegmask)
#' 

bandwidths <- c(1.6,1.9,3.2) #in units of img
lacgb <- function(img,bandwidths){
  if(img$xstep != img$ystep){print("ERROR: image pixels must be square")}
  b <- round(bandwidths/img$xstep)
  b <- unique(bandwidthspixX) 
  bandwidths <- b*img$xstep

  lacs <- mapply(lacgb0,bX=b,bY=b,b=bandwidths,MoreArgs=list(img=img),SIMPLIFY=FALSE)

  return(lacs)
}

lacgb0(img,bX,bY,bX*0.8) #b is bandwidth in img units for the RS correction
lacgb0 <- function(img,bX,bY,b){
  mat <- matrix(1/((1+2*bX)*(1+2*bY)),ncol=1+2*bX,nrow=1+2*bY)
  kernelfcn <- im(mat,xcol=(-bX:bX)*img$xstep,yrow=(-bY:bY)*img$ystep)
 
  areafracs <- convolve.im(img,kernelfcn)

#lacunarity if the box centeres can be everywhere (aka no boundary correction)
  smA <- mean(areafracs) #sample mean
  ss2A <- mean(areafracs^2) #biased sample second moment
  lacA <- ss2A/(smA^2) -1
  if (is.empty(erosion(Window(img),b))){return(list(lacA=lacA,lacRS=NULL))}
  areafracsRS <-  as.im(convolvedimg,W=erosion(Window(img),b))
  smRS <- mean(areafracsRS) #sample mean
  ss2RS <- mean(areafracsRS^2) #biased sample second moment
  lacRS <- ss2RS/(smRS^2) -1
  return(list(lacA=lacA,lacRS=lacRS))
}
lacA
lacRS
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

#TO DO:
## catch warnings about empty RS window and print something more understandble
