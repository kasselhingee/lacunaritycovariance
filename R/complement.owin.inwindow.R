
#' @title The complement of a set within another set
#' @export complement.owin.inwindow
#' @description 
#' Primarily for estinating the complement of a set within an observation window.
#' Given a set \eqn{X} and window \eqn{W} it returns
#' \deqn{W \cap X^c}
#' where $X^c$ is the complement of $X$.



#' @param x A set in owin format (can be either raster or polygonal).
#' @param w An outer window which will contain the complement. Typically an observation window when x is an observed set.
#' @return An owin object.


#' @examples
#' xi <- heather$coarse
#' xic <- complement.owin.inwindow(xi, Frame(heather$coarse))

#tests
# library(testthat)
# w <- disc(radius = 5, centre = c(5,10))
# xi <- heather$coarse
# xic <- complement.owin.inwindow(xi, w)
# expect_equal(area.owin(union.owin(intersect.owin(xi, w), xic)), area.owin(w), tolerance = area.owin(w)/1000)
# expect_equal(area.owin(intersect.owin(xi, xic)), 0)
# xip <- as.polygonal(xi)
# xipc <- complement.owin.inwindow(xip, w)
# expect_equal(area.owin(union.owin(intersect.owin(xip, w), xipc)), area.owin(w), tolerance = area.owin(w)/1000)
# expect_equal(area.owin(intersect.owin(xip, xipc)), 0, tolerance = area.owin(w)/1E6)
# xic <- complement.owin.inwindow(square(), square(r=0.5))

#' @details 
#' The usual complement.owin function in spatstat does not (to my knowledge) accomplish the task that complement.owin.inwindow does as
#' the spatstat function requires a rectangular outer window.

complement.owin.inwindow <- function( x, w){
  xc <- setminus.owin(w, x)
  return(xc)
}
