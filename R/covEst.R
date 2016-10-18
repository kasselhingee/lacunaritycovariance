#' @title Covariance estimates, also known as `two-point probabilities' for stationary RACS
#' @export covariance covarianceEstAtPoint covarianceMapEst_direct
#' @description 
#' These functions estimate the covariance of a stationary random closed set. 
#' The covariance is also known as the two-point coverage probability, and very closely related to the semivariogram.
#'  The covariance of a vector \eqn{v} is the probability of two points separated by a \eqn{v} being covered by \eqn{\Xi}.
#' @author{Kassel Hingee}



#' @param Xi An observation (in spatstat owin format or for \code{covariance} an \code{im} object is allowed) of the RACS of interest.
#' @param w The observation window in \code{owin} format. If itsn't included its taken to be the smallest rectangle enclosing \code{Xi}.
#' @param w The boundary of the observation in OWIN format. This is needed for the cases where the observation is not a rectangular region.
#' @param setCovBoundaryThresh to avoid instabilities of dividing by very small areas, any vector \eqn{v} set covariance of the boundary smaller than this threshold is given a covariance of NA 
#' @param inclraw If TRUE the output will be two \code{im} objects one for the standard estimator and one raw estimate.

#' @return 
#' \item{comp1 }{An estimate (assuming stationarity of \eqn{\Xi}) that two points separated by \eqn{v} will be in \eqn{\Xi}.}
#' \item{comp2 }{Denominator - The set covariance of the boundary \code{w}}
#' \item{comp3 }{Numerator - The set covariance of \eqn{\Xi_{obs}}}
#' @return \code{SpatStat} \code{im} objects containing the estimated covariances.
#' 
#' For \code{covarianceEstAtPoint} these are single numerical values; for \code{covarianceMapEst_direct} they are matrices; for \code{covariance} they are objects of SpatStat's \code{im} class.


#' @examples
#' XiOWIN <- heather$coarse
#' windowOWIN <- Frame(heather$coarse)
#' 
#' data(balcattapark_coarse)
#' xi <- balcattapark_coarse$vegmask
#' w <- balcattapark_coarse$boundary
#' covar <- covariance(xi,w,inclraw=FALSE)
#' plot(covar,clipwin=owin(xrange=c(-10,10),yrange=c(-10,10)),axes=TRUE)
#' 
#' coverageProb <- coveragefrac(XiOWIN,windowOWIN)
#' 
#' covariancePt <- covarianceEstAtPoint(XiOWIN,windowOWIN,c(3,5))
#' 
#' covarianceDirectEst <- covarianceMapEst_direct(XiOWIN,windowOWIN,1,1)
#' filled.contour(covarianceDirectEst$covariance)
#' 
#' 
#' covarianceFcn <- covariance(XiOWIN,windowOWIN)
#' plot(covarianceFcn$covariance)
#' plot(covarianceFcn$covariance - coverageProb*coverageProb)

#' @note 
#' The name of this function is probably a bit confusing. Perhaps call it two-point coverage probability, and the other functions `two-point coverage probability' functions.
#' 
#' Double check that I'm calculating set covariance directly properly (I don't think there is a need to reflect the set, but maybe I got things wrong)
#' 



#' @keywords spatial nonparametric

#' @details 
#' The reduced sample estimator is [1]
#' \eqn{ \hat{C}(v) = \frac{|\Xi \cap W \cap ((\Xi \cap W ) \oplus v)|}{|W \cap W\oplus v|}}
#' and the raw estimate (if requested) is
#' \eqn{\hat{\tilde{C}}(v) = \frac{|\Xi \cap W \cap ((\Xi \cap W ) \oplus v)|}{|W|}.}
#' \code{covariance} uses Fourier transforms to calculate set covariances (using \code{\link[spatstat]{setcov}} function). It is much faster (500 times faster in one comparison) than \code{covarianceMapEst_direct}.
#' Vectors with small set covariance of the window are eliminated (using \code{setCovBoundaryThresh} because they cause the covariance to be enormous)

#' @references [1] Molchanov, I. (1997) Statistics of the Boolean Model for Practitioners and Mathematicians. John Wiley & Sons.

covariance <- function(Xi,w=NULL,inclraw=FALSE,setCovBoundaryThresh = 0.1*area.owin(w)){
   if (is.owin(xi)){
      if (!is.null(w)){xi <- intersect.owin(xi,w)}
      else {w <- Frame(xi)}
      setcovXi <- setcov(xi)
      setcovB <- setcov(w)
   }
   else if (is.im(xi)){
      if (!is.null(w)) {
          winim <- as.im(w,xy=xi)
          xi <- eval.im(xi*winim)
      }
      else {w <- Frame(xi)}
      setcovXi <- imcov(xi)
      setcovB <- setcov(w)
   }
   else {
      print("ERROR: Input xi is not an image or owin object")
      return(NULL)
   }
   setcovB[setcovB<setCovBoundaryThresh] <- NA #make NA any values that are too small and lead to division to close to 0
   covar <- eval.im(setcovXi/setcovB,harmonize=TRUE)
   if (!inclraw) {return(covar)}
   if (inclraw) {
      return(list(rs = covar,
                  raw = setcovXi/area(w)))
   }
}

#catch the warning message about harmonising 

########################################
#' @rdname covariance
#' @details \code{covarianceEstAtPoint} estimates the covariance of a single vector \eqn{v} by ratioing the set covariance of \code{Xi} to the set covariance of of observation window \code{w}. Set covariance is calculated by intersecting a set with a translated copy of itself.
#' @param v A 2D vector in c(x,y) format.

covarianceEstAtPoint <- function(Xi,w,v){
  denominator <- area.owin(intersect.owin(w,shift.owin(w,vec=v)))#need to handle denominator of 0
  if (denominator == 0){
    covarianceEst <- NA
    return(list(numerator = NA, denominator = NA, covarianceEst = NA))
    }
  else {  
    Xiinside <- intersect.owin(Xi,w) 
    numerator <- area.owin(intersect.owin(Xiinside,shift.owin(Xiinside,vec=v)))
    covarianceEst <- numerator/denominator
    }
  
  return(list(numerator = numerator, denominator = denominator, covarianceEst = covarianceEst))
}


#######################################
#' @rdname covariance
#' @details 
#' \code{covarianceMapEst_direct} estimates covariance on a regular grid using the resolution of \code{Xi}. The regular grid extends to \code{maxXshiftdistance} and \code{maxYshiftdistance} in the x and y components respectively. It uses \code{covarianceEstAtPoint} to estimate the covariance at each grid point.
#' Ignores point estimates that use an area smaller than 10% of the window.
#' 
#' @param maxXshiftdistance the maximum size of x-component of vectors \eqn{v} to estimate
#' @param maxYshiftdistance the maximum size of y-component of vectors \eqn{v} to estimate
covarianceMapEst_direct <- function(Xi,w,maxXshiftdistance,maxYshiftdistance){
  windowArea <- area.owin(w)
  #create the vectors for testing
  shiftVectorX <- c(-rev(seq(0,maxXshiftdistance,by=Xi$xstep)),seq(Xi$xstep,maxXshiftdistance,by=Xi$xstep))
  shiftVectorY <- c(-rev(seq(0,maxYshiftdistance,by=Xi$ystep)),seq(Xi$ystep,maxYshiftdistance,by=Xi$ystep))
  covarMap = matrix(nrow=length(shiftVectorY),ncol=length(shiftVectorX))
  numeratorMap = matrix(nrow=length(shiftVectorY),ncol=length(shiftVectorX))
  denominatorMap = matrix(nrow=length(shiftVectorY),ncol=length(shiftVectorX))
  for (i in 1:length(shiftVectorX)){
    for (j in 1:length(shiftVectorY)){
      covarianceEstimation <- covarianceEstAtPoint(Xi,w,c(shiftVectorX[i],shiftVectorY[j]))
      if (is.na(covarianceEstimation$denominator) || (covarianceEstimation$denominator < 0.1*windowArea)){
        covarMap[i,j] <- NA
      }
      else {covarMap[i,j] <- covarianceEstimation$covarianceEst}
      numeratorMap[i,j] <- covarianceEstimation$numerator
      denominatorMap[i,j] <- covarianceEstimation$denominator
    }
  } 
  
  covarianceMap = list(covariance = covarMap, numerator=numeratorMap, denominator = denominatorMap, xcol = shiftVectorX, ycol = shiftVectorY, xstep=Xi$xstep,ystep=Xi$ystep)
  return(covarianceMap) 
}


