#' @title Variance Estimates for Observed Area - v2d
#' @export sae.v2d.mean sae.v2d.var  sae.v2d.wsu.mean  sae.v2d.wsu.var
#' @description Estimates the variance of the area of a cover type observed in a thematic map created using a fallible classifier from remote sensing.
#' @author{Kassel Hingee}



#' @param xi.area Observed area of the class of interest.
#' @param window.area Area of the observation window (do not include NA pixels)
#' @param delta Width of a pixel
#' @param p21 Probability that a location fallibly classified into the cover of interest is really outside the cover of interest
#' @param p12 Probability that a location fallibly classified into the non-interesting class is really the class of interest.
#' @param radius A radius of dependence - see details for meaning.

#' @details 
#' The model for the area in the class of interest is
#' \deqn{|X_1| = |\hat{X}_1|- \Delta^2 \pi r^2 \text{Bin}(n_1,p_{21}) + \Delta^2 \pi r^2 \text{Bin}(n_2,p_{12})}
#' where \eqn{\Delta} is the width of a pixel, 
#' \eqn{p_{21}} and \eqn{p_{12}} are the probability of an arbitrary point classified into class 1 or 2 being misclassified respectively (class 1 is the class of interest, class 2 is everything else),
#' and 
#' \eqn{n_1 = \frac{|\hat{X}_1|}{\pi r^2}},
#' \eqn{n_2 = \frac{|\hat{X}_2|}{\pi r^2}}
#' where \eqn{r} is a guess of the radius of dependence in the error process (e.g.~Pete chose \eqn{5m}) for \eqn{r}.
#' This is like assuming the pixels are the size of \eqn{\pi r^2}, or that error process is random on a version of the image partitioned in cells of area \eqn{\pi r^2}.

#perfect confusion matrix, independent pixels error
sae.v2d.mean <- function(xi.area, window.area, delta, p21, p12, radius){
  return(
    xi.area 
    - ( xi.area * p21 ) 
	+ ( (window.area - xi.area) * p12 )
	)
}

#' @describeIn sae.v2d.mean Expected cover area when including sampling uncertainty
sae.v2d.wsu.mean <- function(xi.area, window.area, delta,  n11, n21, n12, n22, radius){
  return(xi.area * ahat11(n11, n21) + ( (window.area - xi.area) * ahat12(n12, n22)))
}

#' @describeIn sae.v2d.mean  The variance estimate assuming confusion matrix, independent error on undefined sub regions
sae.v2d.var <- function(xi.area, window.area, delta, p21, p12, radius){
  return( delta^2 * (pi * radius^2) * p21 * (1 - p21)
         + delta^2 * (pi * radius^2) * (window.area - xi.area) * p12 * (1 - p12) )
}

#' @describeIn sae.v1b.mean  The variance estimate assuming indepedendent sampling variation in conditional confusion matrix estimate and independent error on undefined subregions
sae.v2d.wsu.var <- function(xi.area, window.area, delta,  n11, n21, n12, n22, radius){
  a11 <- ahat11(n11, n21)
  a12 <- ahat12(n12, n22)
  expectofvar <-  delta^2 * (pi * radius^2) * (1 - a11) * a11
          + delta^2 * (pi * radius^2) * (window.area - xi.area) * (1 - a11) * a11
  varofexpect <- xi.area^2 * ssq11(n11, n21) + (window.area - xi.area)^2 * ssq12(n12, n22)
  return(expectofvar + varofexpect)
}
