#' @title Gliding Box lacunarity from a black and white image
#' @export mvlgb lacgb
#' @importFrom utils installed.packages
#'
#' @description Calculates the gliding box lacunarity
#' @details Calculates the gliding box lacunarity for a given range of box sizes (`radius').
#' The centres are given the pixels of \code{img}. If the package \code{raster} is available a moving window algorithm is may be used (via \code{focal()}), however for now the convolution method in spatstat appears to be faster and is default.
#' 
#' Note: (1) The sidelengths are rounded such that they are an odd number of pixels across. (2) The reduced sample points are given by Minkowski subtraction
#' @references The gliding box algorithm is described in Allain, C. and Cloitre, M. (1991) Characterizing the lacunarity of random and deterministic fractal sets. Physical Review A, 44, 3552-3558.
#'
#' @return An fv object with columns for no edge correction (raw) and reduced sample border correction (RS). The side lengths (labelled \code{s}) are always odd multiples of pixel widths.
#' @param img An image of 0's and 1's, and NA pixels are assumed to be outside the observation window.
#' @param sidelengths A list of suggested box side lengths in the dimensions of \code{img}. Note the actual side lengths used will be the closest multiple of an odd number of pixel widths.
#' @param inclraw If TRUE the function will also return a gliding box lacunarity that ignores edge effects. (Default is FALSE)
#' @param W Optional observation window. The observation window used for the estimator will be the union of \code{W} and the NA pixles in \code{img}.
#' @param method Obsolete.  
#' @examples
#' img <- as.im(heather$coarse,na.replace=0)
#' sidelengths <- c(0.2,1,2.2,10) #in units of img
#' lac <- mvlgb(img,sidelengths, inclraw=TRUE)
#' plot(lac, cbind(RS,raw) ~ s)
#'
mvlgb <- function(img,sidelengths,inclraw=FALSE,W=Frame(img), method=""){
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

  if (!("RcppRoll" %in% installed.packages()[,1])){
     stop("Need RcppRoll installed to calculate gliding box lacunarity")
  }
  
  img[(complement.owin(intersect.owin(W,Frame(img)),frame=Frame(img)))] <- NA  #make sure the pixels outside W are set to NA so that reduce sampling happens naturally ##NOTE: this a time consuming operation that may never be needed
	lacs <- mapply(lacgb0.rcpproll, sidep=2*rpix+1,MoreArgs=list(img=img, W=obsvd),SIMPLIFY=FALSE,inclraw)	
  
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
##The following function calculates lacunarity for a box with sidelengths 2*bX+1 and 2*bY+1 (in pixels). The RS version is automatically calculated by ignoring those boxes that have sums that includa NA values. 
#eg lacgb0(img,5,5,5*0.8)
#the W is only for the raw version and must be an owin object. 
#uses rcpproll
lacgb0.rcpproll <- function(img,sidep,inclraw,W=Frame(img)){
  	mat <- as.matrix(img)
	movline.overrows <- RcppRoll::roll_sum(mat, sidep)
	movline.overrowthencols <- RcppRoll::roll_sum(t(movline.overrows),sidep)*img$xstep*img$ystep
	smRS <- mean(movline.overrowthencols, na.rm=TRUE) #sample mean
	ss2RS <- mean(movline.overrowthencols^2, na.rm=TRUE) #biased sample second moment
	lacRS <- ss2RS/(smRS^2) -1
	
	if (inclraw){
		#lacunarity if the box centeres can be everywhere (aka no boundary correction)
		matPAD <- matrix(0,ncol=ncol(mat)+2*sidep, nrow=nrow(mat)+2*sidep)
		matPAD[1*sidep+1:nrow(img),1*sidep+1:ncol(img)] <- mat
		matPAD[is.na(matPAD)] <- 0
		movline.overrows <- RcppRoll::roll_sum(matPAD, sidep)
		movline.overrowthencols <- RcppRoll::roll_sum(t(movline.overrows),sidep)
		

		numpixelsinwindow <- area.owin(W)/(img$xstep*img$ystep)
		smA <- sum(movline.overrowthencols)*img$xstep*img$ystep/numpixelsinwindow #note this isn't simply the mean of areafracs because the areafracs image is larger than the input image
		ss2A <- sum(movline.overrowthencols^2)*(img$xstep*img$ystep)^2/numpixelsinwindow
		lacA <- ss2A/(smA^2) -1
	}
  if (inclraw){ return(list(raw=lacA,RS=lacRS))}
  else {return(RS=lacRS)}
}

#' @rdname mvlgb
lacgb <- mvlgb

#TO DO:
## catch warnings about empty RS window and print something more understandble
## get lidar data at some point
