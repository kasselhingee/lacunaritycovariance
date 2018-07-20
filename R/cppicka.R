#' @title Picka's Modified Estimator of Coverage Probability
#' @author Kassel Liam Hingee 
#' @import spatstat
#' @export cppicka
#' 
#' @description 
#' Picka2000va modifies coverage probability and covariance estimators to develop estimators of centered covariance and pair-correlation function that are balanced.
#' The entity computed by this function corresponds to Picka's \eqn{\hat{m}_1[\check{0},v]}.
#' The result should be 
#' \deqn{\frac{|\Xi \cap W \cap (W \oplus v)|} {|W \cap (W\oplus v)|} }
#' @param xi An observation of a RACS in \pkg{spatstat}'s \code{owin} or \code{im} format.
#' @param obswin The window of observation (not necessarily rectangular) also in \code{owin} format.
#' @param setcov_boundarythresh Any vector \eqn{v} such that set covariance of the observation window is smaller than this threshold is given a value of NA to avoid instabilities caused by dividing by very small areas 


#' @return An estimate of the coverage probability
#' @details
#' If \code{xi} is in \code{im} format then \code{xi} must be an image of 1s, 0s and NAs
#'  representing inside the set, outside the set and outside the observation window respectively.
#'  \code{coverageprob} will not accept a \code{obswin} argument if \code{xi} is in \code{im} format.
#' 
#' @examples
#' xi <- heather$coarse
#' obswindow <- Frame(heather$coarse)
#' cp <- coverageprob(xi, obswindow)
#' cpp1 <- cppicka(xi, obswindow)
#' cpp1[as.ppp(c(0,0), W = obswindow)]
#' # plot(cpp1[cpp1 < 2, drop = FALSE])
#' 
#' #test that it is doing what we expect for a single vector
#' xi <- shift.owin(square(r = 1), vec = c(2, 2))
#' win <- square(r = 4)
#' #expect cppicka at v = c(2,2) to return 1/(2*2)
#' #expect cppicka at v = c(-2, -2) to return 0
#' cpp1 <- cppicka(xi, obswin = win)
#' cpp1[ppp(x = c(2, -2), y = c(2, -2), window = Frame(cpp1))]

#' @keywords spatial nonparametric
cppicka <- function(xi, obswin = NULL,
        setcov_boundarythresh = NULL) {
  if (is.im(xi)){
    if(!is.null(obswin)){xi[setminus.owin(Frame(xi), obswin)] <- 0}
    #check that xi is only 1s, 0s and NAs
    uvals <- unique(as.list(as.matrix(xi)))
    # convert to images of xi and obswin without NAs
    if ( !all(  (uvals %in% c(0, 1)) | is.na(uvals))  && 
             !all((uvals %in% c(FALSE, TRUE, NA)) | is.na(uvals)) ) {
        stop("Input xi has values other than 0, 1 or NA")
    } else {
        if(!is.null(obswin)){obswin <- as.im(obswin, W = Frame(obswin), xy = xi) * xi}
        else {obswin <- xi}
        #saving all NAs as 0 and everything else as 1 in obswin
        obswin[!is.na(as.matrix(obswin))] <- 1
        obswin[is.na(as.matrix(obswin))] <- 0
        #saving all NAs as 0 and everything else as 1 in xi
        xi[is.na(as.matrix(xi))] <- 0 #turn all NA's in xi to 0s
    }
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

  #now xi and obswin should be of the same style regardless of input
  if (is.null(setcov_boundarythresh)){
    setcov_boundarythresh <- 0.1 * sum(obswin)*obswin$xstep*obswin$ystep
  } else if (setcov_boundarythresh < setcovW$xstep * setcovW$ystep * 1E-8){
    warning("setcov_boundarythresh is smaller than A*1E-8 where A is the size of a pixel.
            This might be smaller than the precision of the set covariance computations.
            Consider setting setcov_boundarythresh higher.")
  }

  #numerator
  top <- convolve.im(xi, obswin, reflectY = TRUE)  #u in (obswin + v) iff (u - v) in obswin. Thus Y needs to have argument reflected
  unitname(top) <- unitname(xi)

  #denominator
  bot <- convolve.im(obswin, obswin, reflectY = TRUE)
  bot[bot < setcov_boundarythresh] <- NA #to remove small denominators
  unitname(bot) <- unitname(obswin)
  return(eval.im( top / bot))
}
