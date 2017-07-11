#' @title A spatial covariance, also known as `two-point probability', estimator for stationary RACS
#' @export covariance
#' @description 
#' These functions estimate the covariance of a stationary RACS. 
#' The covariance is also known as the two-point coverage probability, and very closely related to the semivariogram.
#'  The covariance of a vector \eqn{v} is the probability of two points separated by a \eqn{v} being covered by the random set \eqn{\Xi}
#' \deqn{C(v) = P(\{x,x+v\}\subseteq \Xi).}
#' @author{Kassel Hingee}



#' @param xi An observation of the RACS of interest. It can be in \pkg{spatstat}'s \code{owin} or \code{im} format. If \code{xi} is in \code{im} format then it is assumed that the pixels will be valued 1 (for foreground), 0 (for background) and NA for unobserved.
#' @param w The observation window in \code{owin} format. If itsn't included and \code{xi} is an \code{owin} object then \code{w} is taken to be the smallest rectangle enclosing \code{xi}. If \code{xi} is a \code{im} object than \code{w} is all the non-NA pixels in \code{xi}.
#' @param setCovBoundaryThresh to avoid instabilities of dividing by very small areas, any vector \eqn{v} set covariance of the boundary smaller than this threshold is given a covariance of NA 
#' @param inclraw If TRUE the output will be two \code{im} objects one for the standard estimator and one raw estimate.

#' @return \code{SpatStat} \code{im} objects containing the estimated covariances.


#' @examples
#' xi <- heather$coarse
#' covar <- covariance(xi,inclraw=FALSE)
#' covar <- covariance(as.im(heather$coarse,na.replace=1))

#' @keywords spatial nonparametric

#' @details 
#' The reduced sample estimator is [1]
#' \eqn{ \hat{C}(v) = \frac{|\Xi \cap W \cap ((\Xi \cap W ) \oplus v)|}{|W \cap W\oplus v|}}
#' and the raw estimate (if requested) is
#' \eqn{\hat{\tilde{C}}(v) = \frac{|\Xi \cap W \cap ((\Xi \cap W ) \oplus v)|}{|W|}.}
#' \code{covariance} uses Fourier transforms to calculate set covariances (using \code{\link[spatstat]{setcov}} function). It is much faster (500 times faster in one comparison) than \code{covarianceMapEst_direct}.
#' Vectors with small set covariance of the window are eliminated (using \code{setCovBoundaryThresh} because they cause the covariance to be enormous)

#' @references [1] Molchanov, I. (1997) Statistics of the Boolean Model for Practitioners and Mathematicians. John Wiley & Sons.

covariance <- function(xi,w=NULL,inclraw=FALSE,setCovBoundaryThresh = 0.1*area.owin(w)){
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
      #check that xi is only 1s, 0s and NAs
      uvals <- unique(as.list(as.matrix(xi)))
      if (sum(is.na(uvals), uvals %in% c(0,1)) != length(uvals)){
          print("ERROR: input xi is an image object but has values which aren't 0, 1, or NA")
          return(NULL)
      }
      else { 
          w <- as.owin(xi) #only the non-NA pixels will be in the window
          xi[is.na(as.matrix(xi))] <- 0 #turn all NA's in xi to 0s
      } 
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
