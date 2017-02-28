#' @title Gliding Box lacunarity from a black and white image
#' @export lacgb 
#'
#' @description Calculates the gliding box lacunarity
#' @details Calculates the gliding box lacunarity for a given range of box sizes (`radius').
#' The centres are given the pixels of \code{img}.
#' 
#' Note: (1) The sidelengths are rounded such that they are an odd number of pixels across. (2) The reduced sample points are given by Minkowski subtraction
#' @references The gliding box algorithm is described in Allain, C. and Cloitre, M. (1991) Characterizing the lacunarity of random and deterministic fractal sets. Physical Review A, 44, 3552-3558.
#'
#' @return An fv object with columns for no edge correction (raw) and reduced sample border correction (RS). The side lengths (labelled \code{s}) are always odd multiples of pixel widths.
#' @param img An image of 0's and 1's, and NA pixels are assumed to be outside the observation window.
#' @param sidelengths A list of suggested box side lengths in the dimensions of \code{img}. Note the actual side lengths used will be the closest multiple of an odd number of pixel widths.
#' @param inclraw If TRUE the function will also return a gliding box lacunarity that ignores edge effects.
#' @param W Optional observation window. The observation window used for the estimator will be the union of \code{W} and the NA pixles in \code{img}.
#' @examples
#' img <- as.im(heather$coarse,na.replace=0)
#' sidelengths <- c(0.2,1,2.2,3) #in units of img
#' lac <- lacgb(img,sidelengths)
#' plot(lac, cbind(RS,raw) ~ s)
#'
lacgb <- function(img,sidelengths,inclraw=TRUE,W=Frame(img)){
  if(abs(img$xstep -img$ystep)>1E-2 * img$xstep){print("ERROR: image pixels must be square")}
#convert sidelengths to odd pixel amounts, taking into account that want a distance to edge
  spix <- 1+round((sidelengths-img$xstep)/(2*img$xstep))*2
  spix <- unique(spix)
  rpix <- (spix-1)/2
  sidel <- spix*img$xstep 

#compute observation mask
  if (class(img)!="im"){print("ERROR: input img must be of class im")}
  obsvd <- img
  obsvd[is.finite(img$v)] <- TRUE
  if (class(W)=="im"){obsvd <- eval.im(W * obsvd)}
  obsvd <- as.owin(obsvd) #owin format needed for use of dilation lacgb0 
  if (class(W)=="owin"){obsvd <- intersect.owin(obsvd,W)}
  img[!is.finite(img$v)] <- 0

  lacs <- mapply(lacgb0,bX=rpix,bY=rpix,MoreArgs=list(img=img, W=obsvd),SIMPLIFY=FALSE,inclraw)


if (inclraw){
  nobord <- unlist(lapply(lacs, `[[`, 1) )
  RS <- unlist(lapply(lacs, `[[`, 2) )
  lacsdf <- data.frame(s = sidel,raw=nobord,RS=RS)
  lacfv <- fv(lacsdf,argu="s",valu="RS",
           ylab = expression(MVL[gb]),
	   unitname=unitname(img),
           labl = c("Sidelength","raw","RS"),
           desc = c("Sidelengths of boxes", "Gliding Box Lacunarity ignoring edge effects", "Gliding Box Lacunarity that only uses boxes entirely within the observation")
           )
}
else {
  RS <- unlist(lapply(lacs, `[[`, 1) )
  lacsdf <- data.frame(s = sidel,RS=RS)
  lacfv <- fv(lacsdf,argu="s",valu="RS",
           ylab = expression(MVL[gb]),
	   unitname=unitname(img),
           labl = c("Sidelength","RS"),
           desc = c("Sidelengths of boxes", "Gliding Box Lacunarity that only uses boxes entirely within the observation")
           )
}
  return(lacfv)
}


##########################
##The following function calculates lacunarity for a box with sidelengths 2*bX+1 and 2*bY+1 (in pixels). It also calculates the RS by eroding by `b' where b is in UNITS OF THE IMAGE.
#eg lacgb0(img,5,5,5*0.8)
#W MUST be an owin object

lacgb0 <- function(img,bX,bY,inclraw,W=Frame(img)){
  stopifnot(class(W)=="owin")#for dilation
  img <- img[W] #make sure the pixels outside W are set to NA
  distfromCentrePtofCentrePix <- bX*img$xstep+0.5*img$xstep
  mat <- matrix(1,ncol=round(1+2*bX),nrow=round(1+2*bY))
  useraster = ("raster" %in% installed.packages()[,1])
  
  if (useraster){#use raster's fast moving window function (called focal)
  		xiras <- raster::raster(img)
		#using raster's moving window stuff
		areas.movwind <- raster::focal(xiras, mat)
		areas.movwind <- raster::as.array(areas.movwind)[,,1]
		areas.movwind <- as.im(areas.movwind[nrow(areas.movwind):1,], W=img)*img$xstep*img$ystep
		smRS <- mean(areas.movwind, na.rm=TRUE) #sample mean
		ss2RS <- mean(areas.movwind^2, na.rm=TRUE) #biased sample second moment
		lacRS <- ss2RS/(smRS^2) -1
  }
  if (!useraster | inclraw){ #convolve is only needed if raw is requested, or if it is not possible to use raster
	##building the kernel fcn
	  kernelfcn <- im(mat/((1+2*bY)*(1+2*bX)*img$xstep*img$ystep),xcol=(-bX:bX)*img$xstep,yrow=(-bY:bY)*img$ystep)#because convolve.im approximates the integral the weight here must be area not just number of pixels
	 
	  areafracs <- convolve.im(img,kernelfcn) #this map includes all the box centres that intersect img (aka it is bigger than img) - means no buffer is required for raw lacunarity
  }
  if (!useraster){
	  allowedBoxCentres <- erosionAny(W,Frame(kernelfcn))
	  if (is.empty(allowedBoxCentres)){lacRS=NA}
	  else {
		areafracsRS <-  as.im(areafracs,W=Frame(img)) 
		rsW <- as.im(allowedBoxCentres,xy=areafracsRS) 
		areafracsRS <- eval.im(areafracsRS * rsW)
		smRS <- mean(areafracsRS, na.rm=TRUE) #sample mean
		ss2RS <- mean(areafracsRS^2, na.rm=TRUE) #biased sample second moment
		lacRS <- ss2RS/(smRS^2) -1
	  }
  }

  if (inclraw) {
    #lacunarity if the box centeres can be everywhere (aka no boundary correction)
    numpixinimage <- dim(img)[1]*dim(img)[2]
    smA <- sum(areafracs)/numpixinimage #note this isn't simply the mean of areafracs because the areafracs image is larger than the input image
    ss2A <- sum(areafracs^2)/numpixinimage
    lacA <- ss2A/(smA^2) -1
  }
  if (inclraw){ return(list(raw=lacA,RS=lacRS))}
  else {return(RS=lacRS)}
}
  

#TO DO:
## catch warnings about empty RS window and print something more understandble
## get lidar data at some point
