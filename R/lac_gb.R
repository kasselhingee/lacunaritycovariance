#' @title Gliding Box lacunarity from a black and white image
#' @export lacgb 
#'
#' @description Calculates the gliding box lacunarity
#' @details Calculates the gliding box lacunarity for a given range of box sizes (`radius').
#' The centres are given the pixels of \code{img}.
#' @references **Algorithm described in Allain1991ch**
#'
#' @return An fv object with items for no edge correction and reduced sample border correction
#' @param img An image of 0's and 1's.
#' @param bandwidths A list of bandwidths (half the box sidelengths) in the dimensions of \code{img}
#' @examples
#' data(balcattapark_coarse)
#' img <- as.im(balcattapark_coarse$vegmask)
#' bandwidths <- c(1.6,1.9,3.2,5*0.8) #in units of img
#' lac <- lacgb(img,bandwidths)
#' plot(lac, cbind(RS,nobord) ~ b)
#'
lacgb <- function(img,bandwidths){
  if(abs(img$xstep -img$ystep)>1E-6 * img$xstep){print("ERROR: image pixels must be square")}
#convert bandwiths to pixel amounts
  b <- round(bandwidths/img$xstep)
  b <- unique(b) 
  bandwidths <- b*img$xstep

  lacs <- mapply(lacgb0,bX=b,bY=b,b=bandwidths,MoreArgs=list(img=img),SIMPLIFY=FALSE)
  nobord <- unlist(lapply(lacs, `[[`, 1) )
nobord
  RS <- unlist(lapply(lacs, `[[`, 2) )
  lacsdf <- data.frame(b = bandwidths,nobord=nobord,RS=RS)
  lacfv <- fv(lacsdf,argu="b",valu="RS",
           ylab = "lacunarity",
	   unitname=unitname(img))
  return(lacfv)
}


##########################
##The following function calculates lacunarity for a box with sidelengths 2*bX+1 and 2*bY+1 (in pixels). It also calculates the RS by eroding by `b' where b is in UNITS OF THE IMAGE.
#eg lacgb0(img,5,5,5*0.8)

lacgb0 <- function(img,bX,bY,b){
##building the kernel fcn
  mat <- matrix(1/((1+2*bX)*(1+2*bY)),ncol=1+2*bX,nrow=1+2*bY)
  kernelfcn <- im(mat,xcol=(-bX:bX)*img$xstep,yrow=(-bY:bY)*img$ystep)
 
  areafracs <- convolve.im(img,kernelfcn)

#lacunarity if the box centeres can be everywhere (aka no boundary correction)
  smA <- mean(areafracs) #sample mean
  ss2A <- mean(areafracs^2) #biased sample second moment
  lacA <- ss2A/(smA^2) -1
  if (is.empty(erosion(Window(img),b))){return(list(lacA=lacA,lacRS=NULL))}
  areafracsRS <-  as.im(areafracs,W=erosion(Window(img),b))
  smRS <- mean(areafracsRS) #sample mean
  ss2RS <- mean(areafracsRS^2) #biased sample second moment
  lacRS <- ss2RS/(smRS^2) -1
  return(list(nobord=lacA,RS=lacRS))
}
  

#TO DO:
## catch warnings about empty RS window and print something more understandble
## get lidar data at some point
## a test to make sure the function stays working
