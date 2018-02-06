#' @title Variance Estimates for Observed Area - v2d
#' @export sae.v2d.mean sae.v2d.var
#' @description Estimates the variance of the area of a cover type observed in a thematic map created using a fallible classifier from remote sensing.
#' @author{Kassel Hingee}



#' @param xi.area Observed area of the class of interest.
#' @param window.area Area of the observation window (do not include NA pixels)
#' @param delta Width of a pixel
#' @param prob2gvnobs1 Probability that a location fallibly classified into the cover of interest is really outside the cover of interest
#' @param prob1gvnobs2 Probability that a location fallibly classified into the non-interesting class is really the class of interest.
#' @param radius A radius of dependence - see details for meaning.

#' @details 
#' The model for the area in the class of interest is
#' \deqn{|X_1| = |\hat{X}_1|-\Delta^2$\pi r^2$\text{Bin}(n_1,p_{21}) + \Delta^2$\pi r^2$\text{Bin}(n_2,p_{12})}
#' where $\Delta$ is the width of a pixel, 
#' \eqn{p_{21}} and \eqn{p_{12}} are the probability of an arbitrary point classified into class 1 or 2 being misclassified respectively (class 1 is the class of interest, class 2 is everything else),
#' and 
#' \eqn{n_1 = \frac{|\hat{X}_1|}{\pi r^2}},
#' \eqn{n_2 = \frac{|\hat{X}_2|}{\pi r^2}}
#' where $r$ is a guess of the radius of dependence in the error process (e.g.~Pete chose $5m$) for $r$.
#' This is like assuming the pixels are the size of $\pi r^2$, or that error process is random on a version of the image partitioned in cells of area $\pi r^2$.

#perfect confusion matrix, independent pixels error
sae.v2d.mean <- function(xi.area, window.area, delta, prob2gvnobs1, prob1gvnobs2, radius){
  return(
    xi.area 
    - ( xi.area * prob2gvnobs1 ) 
	+ ( (window.area - xi.area) * prob1gvnobs2 )
	)
}

#' @describeIn sae.v2d.mean  The variance estimate assuming confusion matrix, independent pixels error
sae.v2d.var <- function(xi.area, window.area, delta, prob2gvnobs1, prob1gvnobs2, radius){
  return( delta^2 * (pi * radius^2) * prob2gvnobs1 * (1 - prob2gvnobs1)
         + delta^2 * (pi * radius^2) * (window.area - xi.area) * prob1gvnobs2 * (1 - prob1gvnobs2) )
}