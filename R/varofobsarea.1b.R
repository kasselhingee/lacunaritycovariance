#' @title Variance Estimates for Observed Area - v1b
#' @export sae.v1b.mean sae.v1b.var
#' @description Estimates the variance of the area of a cover type observed in a thematic map created using a fallible classifier from remote sensing.
#' @author{Kassel Hingee}



#' @param xi.area Observed area of the class of interest.
#' @param window.area Area of the observation window (do not include NA pixels)
#' @param delta Width of a pixel
#' @param prob2gvnobs1 Probability that a location fallibly classified into the cover of interest is really outside the cover of interest
#' @param prob1gvnobs2 Probability that a location fallibly classified into the non-interesting class is really the class of interest.

#' @details 
#' The model for the area in the class of interest is
#' \deqn{|X_1| = |\hat{X}_1|-\Delta^2\text{Bin}(n_1,p_{21}) + \Delta^2\text{Bin}(n_2,p_{12})}
#' where $\Delta$ is the width of a pixel, $n_1$ and $n_2$ are the number of pixels in $\hat{X}_1$ and $\hat{X}_2$ respectively, and
#' \eqn{p_{21}} and \eqn{p_{12}} are the probability of an arbitrary point classified into class 1 or 2 being misclassified respectively (class 1 is the class of interest, class 2 is everything else).



#perfect confusion matrix, independent pixels error
sae.v1b.mean <- function(xi.area, window.area, delta, prob2gvnobs1, prob1gvnobs2){
  return(xi.area - ( xi.area * prob2gvnobs1 ) + ( (window.area - xi.area) * prob1gvnobs2))
}

#' @describeIn sae.v1b.mean  The variance estimate assuming confusion matrix, independent pixels error
sae.v1b.var <- function(xi.area, window.area, delta, prob2gvnobs1, prob1gvnobs2){
  return( delta^2 * prob2gvnobs1 * (1 - prob2gvnobs1) + delta^2 * (window.area - xi.area) * prob1gvnobs2 * (1 - prob1gvnobs2) )
}