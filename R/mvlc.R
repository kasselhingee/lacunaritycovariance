#' @title Covariance-based calculations of mass variance lacunarity
#' @export lac
#' @export mvl
#' @export mvlc
#'
#' @description Estimates the mass variance lacunarity (MVL) of a stationary RACS from a bi-tonal image, or calculates the MVL from provided covariance (two-point probability) and coverage fraction.

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
#' sidelengths <- seq(0.3,3,by=0.2)
#' plot(lac(sidelengths,covar,p))
#' otherboxes <- list(owin(xrange=c(-0.15,0.15), yrange=c(-0.15,0.15)),
#'                    square(0.3),square(0.5),disc(0.8),square(1))
#' otherboxes.mvl <- mvlc(otherboxes,covar,p)
#' points(c(0.3,0.5,0.8,1),otherboxes.mvl)
#' 
#' #Test on a Boolean Model
#' lambda <- 2.2064E-3
#' discr <- 10
#' w <- owin(xrange=c(0,100),yrange=c(0,100))
#' xi <- rBooleanDetermDiscs(lambda,discr,w)
#' plot(xi)
#' xiimg <- as.im(xi, W=w, eps=c(0.1,0.1), na.replace=0)
#' #estimated lacunarity
#' mvl.est <- mvlc(c(0.5,0.8,1,2,3,4,5,6),xiim=xiimg)
#' plot(mvl.est)
#' #theoretical lacunarity very different because window is small **I think
#' thcovariance <- thcovarDeterministicDiscs(
#'                  xrange=c(-10,10),
#'		    yrange=c(-10,10),
#'		    eps=c(0.1,0.1),lambda,discr)
#' thcoverageprob <- booldetermdiscs_truecoveragefrac(lambda,discr)
#' mvl.th <- mvlc(c(0.5,0.8,1,2,3,4,5,6),thcovariance,thcoverageprob)
#' plot(add=TRUE, mvl.th, col="red", lty="dashed")

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
     boxcov <- lapply(boxes,setcovsquare,xy=covariance) #theoretical 
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
   A <- A[integrationregion]
   B <- B[integrationregion]
   prdimg <- eval.im(A*B,harmonize=TRUE)
   return(sum(prdimg[,],na.rm=na.rm)*prdimg$xstep*prdimg$ystep)
}

#for a square the set covariance can be calculated analytically using sidelengths
setcovsquare <- function(side,xy=NULL){
  if (is.im(xy)){
     xcol <- xy$xcol
     yrow <- xy$yrow
  }
  if (is.null(xy)) {#defaul to 100 points
     xcol <- seq(-side,side,length.out=100)
     yrow <- seq(-side,side,length.out=100)
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
