#' @title Picka's Reduced Window Estimator of Coverage Probability
#' @author Kassel Liam Hingee 
#' @import spatstat
#' @export cppicka
#' 
#' @description 
#' This function provides estimates of coverage probability from subsets of the observation window,
#'  which are a key component of balanced estimators of covariance, centred covariance, pair-correlation and gliding box lacunarity.
#' @param xi An observation of a RACS of interest as a full binary map (as an \code{im} object) or as the foreground set (as an \code{owin} object).
#' In the latter case the observation window, \code{obswin}, must be supplied.
#' @param obswin If \code{xi} is an \code{owin} object then \code{obswin} is an
#'   \code{owin} object that specifies the observation window.
#' @param setcov_boundarythresh To avoid instabilities caused by dividing by very small quantities, if the set covariance of the observation window
#'  is smaller than \code{setcov_boundarythresh}, then the returned pixel value is NA.


#' @return An \code{im} object. Pixel values correspond to estimates of the coverage probability
#' from the subregion of the observation window, \eqn{W}, that is the intersection of \eqn{W} and \eqn{W} shifted by vector \eqn{v}, where \eqn{v} is the pixel location.
#' 
#' @details
#' The plug-in moment covariance estimator (\code{\link{plugincvc}}) uses less of the observation window than the usual coverage probability estimators.
#' Picka (1997, 2000) created new 'balanced' estimators of centred covariance and pair-correlation
#' that accounted for this difference.
#' A key component of Picka's estimators is an estimate of the coverage probability from the subregion of the binary map that is
#' the intersection between \eqn{W} and \eqn{W} shifted by vector \eqn{v}, where \eqn{W} is the observation window (p.~687, Picka, 2000).
#' If we treat \eqn{X} and \eqn{W} as indicator functions representing the foreground and observation window respectively,
#'  this coverage probability estimator used by Picka is
#' \deqn{ \frac{\int X(u) W(u) W(u - v) du} {\int W(u) W(u - v) du}. }{integral(X(u) W(u) W(u - v) du)  /  integral(W(u) W(u - v) du).}
#' 
#' \code{cppicka} produces these estimates for an array of vectors \eqn{v} using fast Fourier transforms.

#' @references
#' Picka, J.D. (1997) \emph{Variance-Reducing Modifications for Estimators of Dependence in Random Sets}. Ph.D.: Illinois, USA: The University of Chicago.
#' 
#' Picka, J.D. (2000) Variance reducing modifications for estimators of standardized moments of random sets. \emph{Advances in Applied Probability}, 32, 682-700.

#' @examples
#' xi <- heather$coarse
#' obswindow <- Frame(heather$coarse)
#' cp <- coverageprob(xi, obswindow)
#' cpp1 <- cppicka(xi, obswindow)

#' @keywords spatial nonparametric
cppicka <- function(xi, obswin = NULL,
        setcov_boundarythresh = NULL) {
  if (is.im(xi)){
    if(!is.null(obswin)){xi[setminus.owin(Frame(xi), obswin)] <- 0}
    #check that xi is only 1s, 0s and NAs
    isbinarymap(xi, requiretrue = TRUE)
    if(!is.null(obswin)){obswin <- as.im(obswin, W = Frame(obswin), xy = xi) * xi}
    else {obswin <- xi}
    #saving all NAs as 0 and everything else as 1 in obswin
    obswin[!is.na(as.matrix(obswin))] <- 1
    obswin[is.na(as.matrix(obswin))] <- 0
    #saving all NAs as 0 and everything else as 1 in xi
    xi[is.na(as.matrix(xi))] <- 0 #turn all NA's in xi to 0s
  }
  
  if (is.owin(xi)){
    stopifnot(is.owin(obswin))
    if (!is.subset.owin(boundingbox(xi), boundingbox(obswin))){
      warning("Bounding box of xi is not inside bounding box of obswin. This can be caused by xi being a pixel mask, and obswin being a polygon that is not quite aligned with the pixels. The Frame(xi) will be enlarged")
    }
    #the following makes sure convolutions generated when the input frame of xi is much smaller than obswin
    Frame(xi) <- boundingbox(obswin, Frame(xi)) #that obswin and xi are needed on right is when the binary grid in xi doesn't quite match with obswin
    xi <- as.im(xi, W = obswin, na.replace = 0)
    obswin <- as.im(obswin, na.replace = 0, xy = xi)
  }

  #numerator
  top <- convolve.im(xi, obswin, reflectY = TRUE)  #u in (obswin + v) iff (u - v) in obswin. Thus Y needs to have argument reflected
  unitname(top) <- unitname(xi)

  #denominator
  bot <- convolve.im(obswin, obswin, reflectY = TRUE)
  #now xi and obswin should be of the same style regardless of input
  if (is.null(setcov_boundarythresh)){
    setcov_boundarythresh <- 0.1 * sum(obswin)*obswin$xstep*obswin$ystep
  } else if (setcov_boundarythresh < bot$xstep * bot$ystep * 1E-8){
    warning("setcov_boundarythresh is smaller than A*1E-8 where A is the size of a pixel.
            This might be smaller than the precision of the set covariance computations.
            Consider setting setcov_boundarythresh higher.")
  }

  bot[bot < setcov_boundarythresh] <- NA #to remove small denominators
  unitname(bot) <- unitname(obswin)
  return(eval.im( top / bot))
}
