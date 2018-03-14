#' @title Picka's Modified Estimator of Coverage Probability
#' @author Kassel Liam Hingee 
#' @import spatstat
#' @export cppicka
#' 
#' @description 
#' Picka2000va modifies coverage probability and covariance estimators to develop estimators of centered covariance and pair-correlation function that are balanced.
#' The entity computed by this function corresponds to Picka's \eqn{\hat{m}_1[\check{0},v]}.
#' @param xi An observation of a RACS in \pkg{spatstat}'s \code{owin} or \code{im} format.
#' @param obswin The window of observation (not necessarily rectangular) also in \code{owin} format.
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
#' plot(cpp1[cpp1 < 2, drop = FALSE])

#' @keywords spatial nonparametric
cppicka <- function(xi, obswin = NULL){
  if (is.im(xi)){
    stopifnot(is.null(obswin))
    #check that xi is only 1s, 0s and NAs
    uvals <- unique(as.list(as.matrix(xi)))
    # convert to images of xi and obswin without NAs
    if ( !all(  (uvals %in% c(0, 1)) | is.na(uvals))  && 
             !all((uvals %in% c(FALSE, TRUE, NA)) | is.na(uvals)) ) {
        stop("Input xi has values other than 0, 1 or NA")
    } else {
        obswin <- xi
        obswin[as.matrix(xi) == 0] <- 1 #only all 0 or 1 valued pixels will be in obswin
        xi[is.na(as.matrix(xi))] <- 0 #turn all NA's in xi to 0s
    }
  }
  
  if (is.owin(xi)){
    stopifnot(is.owin(obswin))
    xi <- as.im(xi, na.replace = 0)
    obswin <- as.im(obswin, na.replace = 0, xy = xi)
  }

  #now xi and obswin should be of the same style regardless of input
  #numerator
  top <- convolve.im(xi, obswin, reflectY = TRUE)  #u in (obswin + v) iff (u - v) in obswin. Thus Y needs to have argument reflected
  unitname(top) <- unitname(xi)
  bot <- convolve.im(obswin, obswin, reflectY = TRUE)
  unitname(bot) <- unitname(obswin)
  return(eval.im( top / bot))
}
