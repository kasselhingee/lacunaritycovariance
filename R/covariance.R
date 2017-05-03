#' @title Covariance estimates, also known as `two-point probabilities' for stationary RACS
#' @export covariance
#' @description 
#' These functions estimate the covariance of a stationary random closed set. 
#' The covariance is also known as the two-point coverage probability, and very closely related to the semivariogram.
#'  The covariance of a vector \eqn{v} is the probability of two points separated by a \eqn{v} being covered by \eqn{\Xi}.
#' @author{Kassel Hingee}



#' @param xi An observation (in spatstat owin format or black and white \code{im} object) of the RACS of interest.
#' @param w The observation window in \code{owin} format. If itsn't included its taken to be the smallest rectangle enclosing \code{Xi}.
#' @param setCovBoundaryThresh to avoid instabilities of dividing by very small areas, any vector \eqn{v} set covariance of the boundary smaller than this threshold is given a covariance of NA 
#' @param inclraw If TRUE the output will be two \code{im} objects one for the standard estimator and one raw estimate.

#' @return \code{SpatStat} \code{im} objects containing the estimated covariances.


#' @examples
#' xi <- heather$coarse
#' covar <- covariance(xi,inclraw=FALSE)

#' data(balcattapark_coarse)
#' xi <- balcattapark_coarse$vegmask
#' w <- balcattapark_coarse$boundary
#' covar <- covariance(xi,w,inclraw=FALSE)

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
