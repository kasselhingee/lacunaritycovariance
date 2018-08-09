#' @title A spatial covariance, also known as `two-point probability', estimator for stationary RACS
#' @export tradcovarest
#' @description 
#' This function estimate the covariance of a stationary RACS. 
#' The covariance is also known as the two-point coverage probability, and very closely related to the semivariogram.
#'  The covariance of a vector \eqn{v} is the probability of two points separated by a vector \eqn{v} being covered by the random set \eqn{\Xi}
#' \deqn{C(v) = P(\{x,x+v\}\subseteq \Xi).}
#' @author{Kassel Liam Hingee}



#' @param xi An observation of the RACS of interest. It can be in \pkg{spatstat}'s \code{owin} or \code{im} format. If \code{xi} is in \code{im} format then it is assumed that the pixels will be valued 1 (for foreground), 0 (for background) and NA for unobserved.
#' If \code{xi} is in \code{owin} format take care to consider what exactly the observation window is (if none is supplied then it will be assumed that the observation window is the smallest rectangle enclosing \code{xi}).
#' @param obswin The observation window in \code{owin} format. If it isn't included and \code{xi} is an \code{owin} object then \code{obswin} is taken to be the smallest rectangle enclosing \code{xi}. If \code{xi} is a \code{im} object than \code{obswin} is all the non-NA pixels in \code{xi}.
#' @param setcov_boundarythresh Any vector \eqn{v} such that set covariance of the observation window is smaller than this threshold is given a covariance of NA to avoid instabilities caused by dividing by very small areas.
#' If NULL is supplied (default) then 0.1 * area.owin(obswin) is used
#' @param inclraw If TRUE the output will be two \code{im} objects one for the standard estimator and one raw estimate.

#' @return A \pkg{SpatStat} \code{im} object containing the estimated covariance. The grey scale values in this image represent the covariance for an array of vectors. If the raw version is requested then a list of \code{im} objects is returned.


#' @examples
#' xi <- heather$coarse
#' covar <- tradcovarest(xi, inclraw = FALSE)
#' covar <- tradcovarest(as.im(heather$coarse, na.replace = 0))

#' @keywords spatial nonparametric

#' @details 
#' The reduced sample estimator [1] is
#' \deqn{ \hat{C}(v) = \frac{\gamma_{W\cap X}(v)}{\gamma_W(v)}}
#' where \eqn{\gamma_{W}(v)} is the set covariance of the observation window \eqn{|W \cap (W\oplus v)|} and \eqn{\gamma_{W\cap X}(v)} is the set covariance of an observation, \eqn{X}, of the RACS \eqn{\Xi}, 
#' \eqn{\gamma_{W\cap X}(v) = |W \cap X \cap ((W\cap X ) \oplus v)|}.
#' The raw estimate (if requested) is
#' \deqn{\hat{\tilde{C}}(v) = \frac{\gamma_{W\cap X}(v)}{|W|}.}
#' \code{tradcovarest} uses Fourier transforms to calculate set covariances (using the \code{\link{setcov}} function from \pkg{spatstat}). 
#' Vectors with small \eqn{\gamma_W(v)} are eliminated using \code{setcov_boundarythresh} because they cause numerical instabilities.
# I suspect this instabilities are because the fourier transforms are only approximately correct, and 0 is within the approximately correct range.

#' @references [1] Serra, J.P. (1982) Image Analysis and Mathematical Morphology. London; New York: Academic Press.

tradcovarest <- function(xi,
        obswin = NULL,
        inclraw =FALSE,
        setcov_boundarythresh = NULL) {
  if (is.owin(xi)) {
    if (!is.null(obswin)) {xi <- intersect.owin(xi, obswin)}
    else {obswin <- Frame(xi)}
    Frame(xi) <- Frame(obswin)
    setcovxi <- setcov(xi)
    unitname(setcovxi) <- unitname(xi)
    setcovwindow <- setcov(obswin, eps = c(setcovxi$xstep, setcovxi$ystep))
    unitname(setcovwindow) <- unitname(obswin)
  } else if (is.im(xi)) {
    if (!is.null(obswin)) {
        winim <- as.im(obswin, xy = xi)
        xi <- eval.im(xi * winim)
    }
    #check that xi is only 1s, 0s and NAs
    isbinarymap(xi, requiretrue = TRUE)

    obswin <- as.owin(xi) #only the non-NA pixels will be in the window
    xi[is.na(as.matrix(xi))] <- 0 #turn all NA's in xi to 0s
    setcovxi <- imcov(xi)
    unitname(setcovxi) <- unitname(xi)
    setcovwindow <- setcov(obswin, eps = c(setcovxi$xstep, setcovxi$ystep))
    unitname(setcovwindow) <- unitname(obswin)
  }
  else {
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
  if (!inclraw) {return(covar)}
  if (inclraw) {
    return(list(rs = covar,
               raw = setcovxi / area(obswin)))
   }
}

#catch the warning message about harmonising 
