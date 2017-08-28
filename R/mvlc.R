#' @title Covariance-based calculations of mass variance lacunarity
#' @export lac
#' @export mvl
#' @export mvlc
#'
#' @description Estimates the mass variance lacunarity (MVL) of a stationary RACS from a bi-tonal image, or calculates the MVL from provided covariance (two-point probability) and coverage probability. 

#' @details
#' If we denoted the estimated covariance by \eqn{\hat{C}(v)} and coverage probability \eqn{\hat{p}} then the estimate of MVL is
#' \deqn{\frac{1}{\hat{p}^2 |B|^2}\int \gamma_B(v)\hat{C}(v)dv -1 }

#' @param boxes Either a list of sidelengths for square boxes or a list of \code{owin} objects of shape.
#' @param covariance  A \code{im} object containing the covariance function (typically estimated by \code{covariance})
#' @param p The coverage probability. Typically estimated by the fraction of the observation window covered by the set of interest.
#' @param xiim An observation of a stationary RACS in \code{im} format. \code{xiim} must have values of either 1, 0 or NA; 1 denotes inside the RACS, 0 denotes outside, and NA denotes unobserved.

#' @return Either an \code{fv} object containing the MVL values and box side lengths or if \code{boxes} is a list of owin objects then \code{lac} returns a list of corresponding MVL values. Note if NA or NaN values in the \code{covariance} object are used then function will return NA or NaN instead of an MVL value. 

#' @examples
#' xi <- heather$coarse
#' covar <- covariance(xi,inclraw=FALSE)
#' p <- area(xi)/area(Frame(xi))
#' sidelengths <- seq(0.3,14,by=0.2)
#' plot(lac(sidelengths,covar,p))
#' # what is the MVL estimates for boxes that are discs?
#' discboxes <- lapply(sidelengths/2,disc)
#' discmvls <- mvlc(discboxes,covar,p)
#' points(sidelengths,discmvls)
#' 

mvlc <- function(boxes, covariance=NULL, p=NULL, xiim=NULL){
   if (!(is.null(covariance) | is.null(p))){
      if (!is.null(xiim)){cat("WARNING: covariance, p and observation image, xiim, given. Only the covariance and p will be used\n")}
      lacv <- lac.cov(boxes,covariance,p)
      unitname <- unitname(covariance)
   } else if (!is.null(xiim)){
      if (is.null(covariance) & is.null(p) != TRUE){ cat("WARNING: xiim supplied and only one of covariance or p supplied so using xiim\n")}
      p <- sum(xiim)/sum(is.finite(xiim$v))
      w <- as.owin(xiim) #w is observation window - only the non NA values end up in window
      xiim[is.na(xiim$v)] <- 0
      covar <- covariance(xiim,w = w)
      lacv <- lac.cov(boxes, covar, p)
      unitname <- unitname(xiim)
   } else {
      stop("Input requires specification of xiim or covariance and p\n")
   }
 
   if (mode(boxes) %in% c("integer","numeric")){
      lacfv <- fv(data.frame(s = boxes, MVL=lacv),
                  argu="s",
		  valu="MVL",
		  ylab=expression(MVL),
		  unitname=unitname,
		  labl = c("Box Side Length", "MVL"),
		  desc = c("Side length of boxes", "MVL derived from covariance")
		  )
       return(lacfv)
   } else (return(lacv)) 
}

lac.cov <- function(boxes, covariance, p){
  if (mode(boxes) %in% c("integer","numeric")){
     boxcov <- lapply(boxes,setcovsquare) #theoretical 
     boxarea <- boxes^2
  }
  else { #box must be a list of owin objects
     boxcov <- lapply(boxes,setcov) #numerical
     boxarea <- lapply(boxes,area.owin)
     boxarea <- unlist(boxarea)
  }

  integrationresults <- mapply(innerprod.im,boxcov,list(covariance),na.rm=FALSE,SIMPLIFY=FALSE)# the list around the covariance is necessary to stop mapply unlisting the image itself

  lac <- unlist(integrationresults)/(p^2 *boxarea^2) -1
  return(lac)
}


innerprod.im <- function(A,B,na.rm=FALSE){
   integrationregion <- intersect.owin(Frame(A),Frame(B))
   #got to do the harmonisation manually so that NA values that the subsetting operation doesn't introduce NA values
   harmgrid <- as.mask(integrationregion, eps=c(min(A$xstep,B$xstep),min(A$ystep,B$ystep)))
   A2 <- as.im(A, xy=harmgrid)
   B2 <- as.im(B, xy=harmgrid)
   prdimg <- eval.im(A2*B2,harmonize=FALSE)
   return(sum(prdimg[,],na.rm=na.rm)*prdimg$xstep*prdimg$ystep)
}
#tests of innerprod.im:
#innerprod.im(as.im(square(1)),as.im(square(1),value=1))
#natest:
#imna <- as.im(square(1.01),value=1)
#imna[as.ppp(c(0.5,0.5),W=Frame(imna))] <- NA
#innerprod.im(as.im(square(1)),imna, na.rm=FALSE)
#innerprod.im(as.im(square(1)),imna, na.rm=TRUE)

#innerprod.im(as.im(square(1)),as.im(square(0.25),value=1))

#should be close to 0 (orthogonal):
#innerprod.im(as.im(function(x,y) {sin(x)},W=square(7*pi),eps=0.01),as.im(function(x,y) {sin(2*x)},W=square(2*pi),eps=0.01))
#should be very non-zero
#innerprod.im(as.im(function(x,y) {sin(x)},W=square(7*pi),eps=0.01),as.im(function(x,y) {sin(x)},W=square(2*pi),eps=0.01))
#it should be (and is) equal to this: sum(as.im(function(x,y){sin(x)*sin(x)},W=square(2*pi),eps=0.01))*0.01*0.01

#for a square the set covariance can be calculated analytically using sidelengths
setcovsquare <- function(side,xy=NULL){
  if (is.im(xy)){
     xcol <- xy$xcol
     yrow <- xy$yrow
  }
  if (is.null(xy)) {#defaul to 100 points
     step <- side/2^7
     xcol <- seq(-side-step,side+step,by=step)
     yrow <- seq(-side-step,side+step,by=step)
  }
  #caculate set covariance for each of the vectors given by xcol and yrow. It could 1/4 more efficient by only looking at the positive quadrant and reflect IF the vectors of interest are symmetrical
  #for a box with side length s the set covariance of (x,y) is:
  #(s-|x|)*(s-|y|) whilst either side is 0
  xsize <- pmax(side-abs(xcol),0)
  xcol <- xcol[xsize >0] #save the columns that lead to non-zero image before deleting these columns.
  xsize <- xsize[xsize>0]
  ysize <- pmax(side-abs(yrow),0)
  yrow <- yrow[ysize >0] #save the rows that lead to non-zero image before deleting these columns.
  ysize <- ysize[ysize>0]
  boxcovV <- outer(ysize,xsize) #first element becomes the rows.
  return(im(boxcovV,xcol = xcol, yrow =yrow))
}

#' @rdname mvlc 
mvl <- mvlc

#' @rdname mvlc
lac <- mvlc
