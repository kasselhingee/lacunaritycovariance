#' @title Covariance (Two-Point Probability) of a Stationary RACS
#' @export covariance
#' @description Estimates covariance using the reduced sample method suggested by Molchanov ... and a raw uncorrected estimate if requested.

#' @details The reduced sample estimator is [1]
#' \eqn{ \hat{C}(v) = \frac{|\Xi \cap W \cap ((\Xi \cap W ) \oplus v)|}{|W \cap W\oplus v|}}
#' and the raw estimate is
#' \eqn{\hat{\tilde{C}}(v) = \frac{|\Xi \cap W \cap ((\Xi \cap W ) \oplus v)|}{|W|}.}
#' @references [1] Molchanov, I. (1997) Statistics of the Boolean Model for Practitioners and Mathematicians. John Wiley & Sons.

#' @param Xi An observed RACS in an \code{owin} object or black and white image (\code{im} format).
#' @param w The observation window in \code{owin} format. If itsn't included its taken to be the smallest rectangle enclosing \code{Xi}.
#' @param inclraw If TRUE the output will be two \code{im} objects one for the standard estimator and one raw estimate.

#' @return \code{SpatStat} \code{im} objects containing the estimated covariances.

#' @examples
#' data(balcattapark_coarse)
#' xi <- balcattapark_coarse$vegmask
#' w <- balcattapark_coarse$boundary
#' covar <- covariance(xi,w,inclraw=FALSE)
#' plot(covar,clipwin=owin(xrange=c(-10,10),yrange=c(-10,10)),axes=TRUE)

covariance <- function(xi,w=NULL,inclraw=FALSE){
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
   covar <- eval.im(setcovXi/setcovB,harmonize=TRUE)
   if (!inclraw) {return(covar)}
   if (inclraw) {
      return(list(rs = covar,
                  raw = setcovXi/area(w)))
   }
      #imcov for xi an image
}

#catch the warning message about harmonising 
#what to do about division by things close to 0?
