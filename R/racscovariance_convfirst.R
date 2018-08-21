#' @title Balanced spatial covariance estimation, also known as `two-point probability', estimator for stationary RACS.
#' @export byconv_cvchats  cvchats_convolves 
#' @description 
#' This function estimate the covariance of a stationary RACS. 
#' A variety of balanced, partially balanced and classical estimates are available. 
#' These functions differ to \code{\link{racscovariance}} by computing separately with the convolutions used in the numerator and denominator until end.

 
#' @param xi An observation of the RACS of interest. It can be in \pkg{spatstat}'s \code{owin} or \code{im} format. If \code{xi} is in \code{im} format then it is assumed that the pixels will be valued 1 (for foreground), 0 (for background) and NA for unobserved.
#' If \code{xi} is in \code{owin} format take care to consider what exactly the observation window is (if none is supplied then it will be assumed that the observation window is the smallest rectangle enclosing \code{xi}).
#' @param obswin The observation window in \code{owin} format. If it isn't included and \code{xi} is an \code{owin} object then \code{obswin} is taken to be the smallest rectangle enclosing \code{xi}. If \code{xi} is a \code{im} object than \code{obswin} is all the non-NA pixels in \code{xi}.
#' @param setcov_boundarythresh Any vector \eqn{v} such that set covariance of the observation window
#'  is smaller than this threshold is given a covariance of NA to avoid instabilities caused by dividing by very small areas, 
#' @param phat The classical estimate of coverage probability,
#'  which is the observed area in \code{xi} divided by the total area of the observation window.
#'  See \code{coverageprob} for more information.
#' @param xy A raster object that specifies the pixel coordinates of the desired covariance image. In the same vein and as.mask in spatstat. **INCORRECT IT IS PASSED TO OBJECTS BEFORE setcov()
#' @param modifications A list of strings specifying desired modifications or functions to apply to cvchat, cpp1 and phat.
#'  modifications = "all" will select all inbuilt modifications. See details. 
#' @param drop If TRUE and one modification selected then the returned value will be a single \code{im} object and not a list of \code{im} object.

#' @param xixi The convolution of a set representing xi with itself
#' @param winwin The convolution of the observation window with itself
#' @param xiwin The convolution of the set xi with the observation window

#' @examples
#' xi <- heather$coarse
#' obswin <- Frame(xi)
#' balancedcvchats <- byconv_cvchats(xi, obswin = Frame(xi), modifications = "all")

#' xixi <- setcov(xi, xy = xi)
#' winwin <- setcov(obswin, xy = xi)
#' xiwin <- setcov(xi, obswin, xy = xi)

byconv_cvchats <- function(xi, obswin,
        xy = NULL,
        setcov_boundarythresh = NULL,
        modifications = "all",
        drop = FALSE){
  if (is.null(setcov_boundarythresh)){
    setcov_boundarythresh <- 0.1 * area.owin(obswin)
  }
  if (is.null(xy)){
    if(is.mask(xi)){
      xy <- xi
    }  else {
      stop("xy must be supplied")
    }
  }
  xixi <- setcov(xi, xy = xy)
  winwin <- setcov(obswin, xy = xy)
  winwin[winwin < setcov_boundarythresh] <- NA #to remove small denominators
  xiwin <- setcov(xi, obswin, xy = xy)
  xiwin[winwin < setcov_boundarythresh] <- NA #to remove small denominators
  phat <- area.owin(xi) / area.owin(obswin)
  
  cvchats <- cvchats_convolves(xixi, winwin, xiwin, phat, modifications = modifications, drop = drop) 
  return(cvchats)
}


#' @describeIn byconv_cvchats Applies multiple modifications simultaneously from a precomputed convolutions xi*xi, w*w, xi*w and phat
cvchats_convolves <- function(xixi, winwin, xiwin = NULL, phat = NULL, modifications = "all", drop = FALSE){
  harmonised <- harmonise.im(xixi = xixi, winwin = winwin, xiwin = xiwin)
  xixi <- harmonised$xixi
  winwin <- harmonised$winwin
  xiwin <- harmonised$xiwin
  fcns <- list(
         trad = cvchat_trad,
         mattfeldt = cvchat_mattfeldt_add,
         pickaint = cvchat_picka_int,
         pickaH = cvchat_picka_H
  )
  if ((modifications == "all")[[1]]) {modifications <- names(fcns)}
  fcnstouse <- fcns[names(fcns) %in% modifications]
  isfunction <- unlist(lapply(modifications, function(x) "function" %in% class(x)))
  modificationsnotused <- modifications[!( (modifications %in% names(fcns)) | isfunction)]
  
  fcnstouse <- c(fcnstouse, modifications[isfunction]) #add user specified modification
  
  if(length(modificationsnotused) > 0){stop(
    paste("The following modifications are not recognised as existing function names or as a function:", modificationsnotused))}
  balancedcvchats <- lapply(fcnstouse, function(x) do.call(x, args = list(xixi = xixi, winwin = winwin, xiwin = xiwin, phat = phat)))
  if (drop & (length(balancedcvchats) == 1)) {return(balancedcvchats[[1]])}
  else { return(as.imlist(balancedcvchats)) }
}



cvchat_trad <- function(xixi, winwin, xiwin = NULL, phat = NULL){
  return(xixi / winwin)
}

cvchat_mattfeldt_add <- function(xixi, winwin, xiwin = NULL, phat = NULL){
  return((xixi - (0.5 * (xiwin + reflect.im(xiwin)))^2 / winwin ) / winwin  +  phat^2)
}

cvchat_picka_int <- function(xixi, winwin, xiwin = NULL, phat = NULL){
  return((xixi - xiwin * reflect.im(xiwin) / winwin) / winwin + phat ^2 )
}

cvchat_picka_H <- function(xixi, winwin, xiwin, phat){
  return((xixi  - phat * xiwin - phat * reflect.im(xiwin)) / winwin  + 2* phat^2 )
}
