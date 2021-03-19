#' @title Plug-in moment covariance estimator
#' @export plugincvc
#' @description 
#' This function computes the plug-in moment covariance estimate of a stationary RACS from a binary map.
#' For a stationary RACS, \eqn{\Xi}, the covariance 
#' for a vector \eqn{v} is the probability of two points separated by a vector \eqn{v} are covered by \eqn{\Xi}
#' \deqn{C(v) = P(\{x,x+v\}\subseteq \Xi).}{C(v) = P({x, x+ v} in Xi).}
#' @author{Kassel Liam Hingee}



#' @param xi An observation of a RACS of interest as a full binary map (as an \code{im} object) or as the foreground set (as an \code{owin} object).
#' In the latter case the observation window, \code{obswin}, must be supplied.
#' @param obswin If \code{xi} is an \code{owin} object then \code{obswin} is an
#'   \code{owin} object that specifies the observation window.
#' @param setcov_boundarythresh To avoid instabilities caused by dividing by very small quantities, if the set covariance of the observation window
#'  is smaller than \code{setcov_boundarythresh}, then the covariance is given a value of NA. 

#' @return A \pkg{SpatStat} \code{im} object containing the estimated covariance.


#' @examples
#' xi <- as.im(heather$coarse, na.replace = 0)
#' covar <- plugincvc(xi)

#' @keywords spatial nonparametric

#' @details 
#' The plug-in moment covariance estimator is (Serra, 1982)
#' \deqn{ \hat{C}(v) = \frac{\gamma_{W\cap X}(v)}{\gamma_W(v)}}{ C(v) = gammaWX(v) / gammaW(v) }
#' where \eqn{\gamma_{W}(v)}{gammaW(v)} is the set covariance of the observation window \eqn{W} 
#' and \eqn{\gamma_{W\cap X}(v)}{gammaWX(v)} is the set covariance of the foreground within \eqn{W}.

#' \code{plugincvc} uses Fourier transforms to calculate the set covariance (using the \code{\link[spatstat.geom]{setcov}} of the foreground and observation window.
#' Vectors with small \eqn{\gamma_W(v)}{ gammaW(v) } are eliminated using \code{setcov_boundarythresh} 
#' as division by small values is numerically unstable.
#' 
# I suspect this instabilities are because the fourier transforms are only approximately correct, and 0 is within the approximately correct range.

#' @references 
#' Serra, J.P. (1982) \emph{Image Analysis and Mathematical Morphology}. London; New York: Academic Press.

plugincvc <- function(xi,
        obswin = NULL,
        setcov_boundarythresh = NULL) {
  if (is.im(xi)) {
    if (!is.null(obswin)) {
        winim <- as.im(obswin, xy = xi)
        xi <- eval.im(xi * winim)
    }
    #check that xi is only 1s, 0s and NAs
    isbinarymap(xi, requiretrue = TRUE)

    obswin <- as.owin(xi) #only the non-NA pixels will be in the window
    xi[is.na(as.matrix(xi))] <- 0 #turn all NA's in xi to 0s
    setcovxi <- imcov(xi)
    setcovwindow <- setcov(obswin, eps = c(setcovxi$xstep, setcovxi$ystep))
  } else if (is.owin(xi)) {
    stopifnot(is.owin(obswin))
    xi <- intersect.owin(xi, obswin)
    Frame(xi) <- Frame(obswin)  #I think rebound.owin, used in Frame<-.owin removes the unitname
    unitname(xi) <- unitname(obswin) #this line fixes the issue of lost unitname until spatstat is updated
    setcovxi <- setcov(xi)
    setcovwindow <- setcov(obswin, eps = c(setcovxi$xstep, setcovxi$ystep))
  } else {
    stop("Input xi is not an image or owin object")
  }
  #make NA any values that are too small and lead to division to close to 0
  if (is.null(setcov_boundarythresh)){
    setcov_boundarythresh <- 0.1 * area.owin(obswin)
  } else if (setcov_boundarythresh < setcovwindow$xstep * setcovwindow$ystep * 1E-8){
    warning("setcov_boundarythresh is smaller than A*1E-8 where A is the size of a pixel.
            This might be smaller than the precision of the set covariance computations.
            Consider setting setcov_boundarythresh higher.")
  }
  setcovwindow[setcovwindow < setcov_boundarythresh] <- NA
  harmims <- harmonise.im(setcovxi, setcovwindow)
  covar <- harmims[[1]]/harmims[[2]]
  return(covar)
}

#catch the warning message about harmonising 
