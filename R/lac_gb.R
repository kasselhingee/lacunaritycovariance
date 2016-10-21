#' @title Gliding Box lacunarity from a black and white image
#' @export lacgb 
#'
#' @description Calculates the gliding box lacunarity
#' @details Calculates the gliding box lacunarity for a given range of box sizes (`radius').
#' The centres are given the pixels of \code{img}.
#' 
#' Note: (1) The bandwidths are rounded such that box sidelengths are an odd number of pixels across. (2) The reduced sample points are given by erosion of the image by a disc with radius bandwidth+0.5*pixelwidth, which is note quite erosion by a square.
#' @references The gliding box algorithm is described in Allain, C. and Cloitre, M. (1991) Characterizing the lacunarity of random and deterministic fractal sets. Physical Review A, 44, 3552-3558.
#'
#' @return An fv object with columns for no edge correction (raw) and reduced sample border correction (RS)
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
  if(abs(img$xstep -img$ystep)>1E-2 * img$xstep){print("ERROR: image pixels must be square")}
#convert bandwiths to pixel amounts
  b <- round(bandwidths/img$xstep)
  b <- unique(b) 
  bandwidths <- b*img$xstep

  lacs <- mapply(lacgb0,bX=b,bY=b,b=bandwidths,MoreArgs=list(img=img),SIMPLIFY=FALSE)
  nobord <- unlist(lapply(lacs, `[[`, 1) )
  RS <- unlist(lapply(lacs, `[[`, 2) )
  lacsdf <- data.frame(b = bandwidths,raw=nobord,RS=RS)
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
  mat <- matrix(1/((1+2*bY)*(1+2*bX)*img$xstep*img$ystep),ncol=round(1+2*bX),nrow=round(1+2*bY))
  kernelfcn <- im(mat,xcol=(-bX:bX)*img$xstep,yrow=(-bY:bY)*img$ystep)
 
  areafracs <- convolve.im(img,kernelfcn) #this map includes all the box centres that intersect img (aka it is bigger than img) - means no buffer is required for raw lacunarity

#lacunarity if the box centeres can be everywhere (aka no boundary correction)
  smA <- mean(areafracs) #sample mean
  ss2A <- mean(areafracs^2) #biased sample second moment
  lacA <- ss2A/(smA^2) -1
  if (is.empty(erosion(Frame(img),b+0.5*img$xstep))){return(list(lacA=lacA,lacRS=NULL))}
  areafracsRS <-  as.im(areafracs,W=erosion(Frame(img),b+0.5*img$xstep)) #note erosion by distance b is not quite the same as erosion by a square of "radius" b
  smRS <- mean(areafracsRS) #sample mean
  ss2RS <- mean(areafracsRS^2) #biased sample second moment
  lacRS <- ss2RS/(smRS^2) -1
  return(list(raw=lacA,RS=lacRS))
}
  

#TO DO:
## catch warnings about empty RS window and print something more understandble
## get lidar data at some point
