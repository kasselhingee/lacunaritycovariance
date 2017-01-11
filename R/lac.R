#' @title Estimates of Theoretical Lacunarity
#' @export lac
#'
#' @description Estimates the Theoretical Lacunarity of a stationary RACS using estimated covariance (two-point probability) and coverage fraction.

#' @details
#' Denoted the estimated covariance by \eqn{\hat{C}(v)} and coverage probability \eqn{\hat{p}} then the estimate lacunarity is
#' \deqn{\frac{1}{\hat{p}^2 |B|^2}\int \gamma_B(v)\hat{C}(v)dv -1 }

#' @param covariance  A \code{im} object containing the covariance function (typically estimated by \code{covariance})
#' @param p The coverage probability. Typically estimated by the fraction of the observation window covered by the set of interest.
#' @param boxes Either a list of sidelengths for square boxes or a list of \code{owin} objects of shape.

#' @return A list of mass-variance lacunarity estimates corresponding to \code{boxes}. If NA or NaN values in the \code{covariance} object are used then function return NA or NaN instead of a lacunarity estimate. 

#' @examples
#' xi <- heather$coarse
#' covar <- covariance(xi,inclraw=FALSE)
#' p <- area(xi)/area(Frame(xi))
#' sidelengths <- c(0.3,0.5,0.8,1)
#' lac(sidelengths,covar,p)
#' otherboxes <- list(square(0.3),square(0.5),disc(0.8),square(1))
#' lac(otherboxes,covar,p)
#' 
#' #Test on a Boolean Model
#' lambda <- 2.2064E-3
#' discr <- 10
#' w <- owin(xrange=c(0,100),yrange=c(0,100))
#' xi <- rBooleanDetermDiscs(lambda,discr,w)
#' plot(xi)
#' xiimg <- as.mask(xi,eps=c(0.1,0.1))
#' #estimated lacunarity
#' estcov <- covariance(xiimg,w=w)
#' estp <- area(xiimg)/area(w)
#' lac(c(0.5,0.8,1,2,3,4,5,6),estcov,estp)
#' #theoretical lacunarity very different because window is small **I think
#' thcovariance <- thcovarDeterministicDiscs(
#'                  xrange=c(-10,10),
#'		    yrange=c(-10,10),
#'		    eps=c(0.1,0.1),lambda,discr)
#' thcoverageprob <- booldetermdiscs_truecoveragefrac(lambda,discr)
#' lac(c(0.5,0.8,1,2,3,4,5,6),thcovariance,thcoverageprob)
#' 
#' #look at why there is a difference
#' estp
#' thcoverageprob
#' plot(as.solist(list(estcov,thcovariance)),axes=TRUE)

lac <- function(boxes, covariance, p){
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


