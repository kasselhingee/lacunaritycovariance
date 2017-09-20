#' @title Gliding Box estimation of mass variance lacunarity
#' @export mvlgb lacgb
#' @importFrom utils installed.packages
#'
#' @description Calculates the gliding box estimate [1] of mass variance lacunarity from a bi-tonal image.
#' @details Calculates the gliding box estimate [1] of mass variance lacunarity for a given range of square box sizes. 
#' The algorithm uses the pixel locations in \code{img} as an array of box centre locations to compute
#'  the mean and variance of the area in a random box of a given size.
#' Locations where the box is not completely within the observation window are ignored.
#'  
#' WARNING: This function needs the \code{roll_sum} function in \code{RcppRoll} to operate.
#' 
#' Note: The side lengths are rounded such that they are an odd number of pixels across.
#' 
#' @references [1] Allain, C. and Cloitre, M. (1991) Characterizing the lacunarity of random and deterministic fractal sets. Physical Review A, 44, 3552-3558.
#'
#' @return An \code{fv} object with a column called \code{MVL} for the usual gliding box estimate described in [1]. 
#'  The side lengths (labelled \code{s}) are always odd multiples of the pixel width.
#'  
#' Another column, called \code{raw}, is included if \code{inclraw=TRUE}. 
#' This column represents a version of the gliding box algorithm that pretends that the items/cover of interest are fully observed
#'  (i.e. that the observation window is the entire plane).
#' @param img An image of pixels valued either \code{0}, \code{1} or \code{NA}. \code{NA} valued pixels are assumed to be outside the observation window.
#' @param sidelengths A list of suggested box side lengths in the same units as \code{img}. Note the actual side lengths used will be the closest multiple of an odd number of pixel widths.
#' @param inclraw If TRUE the function will also return a gliding box estimate that ignores edge effects. (Default is FALSE)
#' @param W Optional observation window. The observation window used for the estimator will be the intersection of \code{W} and the pixels that are not \code{NA} in \code{img}.
#' @examples
#' img <- as.im(heather$coarse,na.replace=0)
#' sidelengths <- seq(0.2,14,by=0.2) #in units of img
#' lac <- mvlgb(img,sidelengths, inclraw=TRUE)
#' plot(lac, cbind(MVL,raw) ~ s)
#'
#' @keywords spatial nonparametric 
mvlgb <- function(img,sidelengths,inclraw=FALSE,W=Frame(img)){
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
     stop("RcppRoll must be installed to calculate gliding box lacunarity")
  }
  
  img[(complement.owin(intersect.owin(W,Frame(img)),frame=Frame(img)))] <- NA  #make sure the pixels outside W are set to NA so that reduce sampling happens naturally ##NOTE: this a time consuming operation that may never be needed
	lacs <- mapply(lacgb0.rcpproll, sidep=2*rpix+1,MoreArgs=list(img=img, W=obsvd),SIMPLIFY=FALSE,inclraw)	
  
if (inclraw){
  nobord <- unlist(lapply(lacs, `[[`, 1) )
  MVL <- unlist(lapply(lacs, `[[`, 2) )
  lacsdf <- data.frame(s = sidel,raw=nobord,MVL=MVL)
  lacfv <- fv(lacsdf,argu="s",valu="MVL",
           ylab = expression(MVL[gb]),
	   unitname=unitname(img),
           labl = c("side length","raw","MVL"),
           desc = c("side lengths of boxes", "Gliding Box Lacunarity ignoring edge effects", "Gliding Box Lacunarity that only uses boxes entirely within the observation")
           )
}
else {
  MVL <- unlist(lapply(lacs, `[[`, 1) )
  lacsdf <- data.frame(s = sidel, MVL=MVL)
  lacfv <- fv(lacsdf,argu="s",valu="MVL",
           ylab = expression(MVL[gb]),
	   unitname=unitname(img),
           labl = c("side length","MVL"),
           desc = c("side lengths of boxes", "Gliding Box Lacunarity that only uses boxes entirely within the observation")
           )
}
  return(lacfv)
}


##########################
##The following function calculates lacunarity for a box with side lengths 2*bX+1 and 2*bY+1 (in pixels). The RS version is automatically calculated by ignoring those boxes that have sums that includa NA values. 
#eg lacgb0(img,5,5,5*0.8)
#the W is only for the raw version and must be an owin object. 
#uses rcpproll
lacgb0.rcpproll <- function(img,sidep,inclraw,W=Frame(img)){
  	mat <- as.matrix(img)
  if ((sidep > nrow(mat)) | (sidep > ncol(mat))){
    mvlgb.rs <- NA
  }
  else {
  	movline.overrows <- RcppRoll::roll_sum(mat, sidep)
  	movline.overrowthencols <- RcppRoll::roll_sum(t(movline.overrows),sidep)*img$xstep*img$ystep
  	sampmean.rs <- mean(movline.overrowthencols, na.rm=TRUE) #sample mean
  	samp2ndmom.rs <- mean(movline.overrowthencols^2, na.rm=TRUE) #biased sample second moment
  	mvlgb.rs <- samp2ndmom.rs/(sampmean.rs^2) -1
  }
	
	if (inclraw){
		#lacunarity if the box centeres can be everywhere (aka no boundary correction)
		matpad <- matrix(0,ncol=ncol(mat)+2*sidep, nrow=nrow(mat)+2*sidep)
		matpad[1*sidep+1:nrow(img),1*sidep+1:ncol(img)] <- mat
		matpad[is.na(matpad)] <- 0
		movline.overrows <- RcppRoll::roll_sum(matpad, sidep)
		movline.overrowthencols <- RcppRoll::roll_sum(t(movline.overrows),sidep)
		

		numpixelsinwindow <- area.owin(W)/(img$xstep*img$ystep)
		sampmean.raw <- sum(movline.overrowthencols)*img$xstep*img$ystep/numpixelsinwindow #note this isn't simply the mean of areafracs because the areafracs image is larger than the input image
		samp2ndmom.raw <- sum(movline.overrowthencols^2)*(img$xstep*img$ystep)^2/numpixelsinwindow
		mvlgb.raw <- samp2ndmom.raw/(sampmean.raw^2) -1
	}
  if (inclraw) {return(list(raw = mvlgb.raw, RS = mvlgb.rs))}
  else {return(RS = mvlgb.rs)}
}

#' @rdname mvlgb
lacgb <- mvlgb

