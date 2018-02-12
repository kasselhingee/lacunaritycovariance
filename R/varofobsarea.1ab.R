#' @title Variance Estimates for Observed Area - v1ab
#' @export sae.v1ab.mean sae.v1ab.var
#' @description Estimates the variance of the area of a cover type observed in a thematic map created using a fallible classifier from remote sensing. 
#' Sampling error in the confusion matrix, assume region is a good representation of population.
#' @author{Kassel Hingee}



#' @param xi.area Observed area of the class of interest.
#' @param window.area Area of the observation window (do not include NA pixels)
#' @param n11 The number of samples that were truly class 1 and fallibly classified as class 1
#' @param n21 The number of samples that were truly class 2 and fallibly classified as class 1
#' @param n12 The number of samples that were truly class 1 and fallibly classified as class 2
#' @param n22 The number of samples that were truly class 2 and fallibly classified as class 2

#' @details 
#' The model for the area in the class of interest is
#' \deqn{
#' hat{|X_1|} = |\hat{X}_1|\hat{a}_{11} + |\hat{X}_2|\hat{a}_{12}
#' }
#' Where \eqn{a_{ij}} are estimates of the entry of the confusion matrix calculated by
#' \deqn{
#' \hat{a}_{ij} &:= \frac{n_{ij}}{n_{.j}}
#' }
#' where \eqn{n_{ij}} is the number of samples fallibly classed \eqn{j} and truly class \eqn{i}, and 
#' \eqn{n_{.j}} denotes \eqn{\sum_i^g n_{ij}}.
#'
#' This model only accounts for the sampling error of the confusion matrix, and thus is kinda like assuming that the region of interest IS the population
#' \deqn{
#' \var\left[\hat{|X_1|}\right] \approx |\hat{X}_1|^2\frac{1}{n_{.1} - 1} \left( n_{11} (1 - \hat{a}_{11})^2 + (n_{.1} - n_{11})(- \hat{a}_{11})^2 \right)
#' \quad \quad + |\hat{X}_2|^2\frac{1}{n_{.2} - 1} \left( n_{12} (1 - \hat{a}_{12})^2 + (n_{.2} - n_{12})(- \hat{a}_{12})^2 \right)
#' }
#'
#' @section WARNING
#' The confusion matrix for the 2009 data has an estimate of 
#' proportion of non-tree that should be tree. This will
#' not be respresented well by a t-distribution?


# sampling error in the confusion matrix.
sae.v1ab.mean <- function(xi.area, window.area, n11, n21, n12, n22){
  return(xi.area * (n11 / (n11 + n21)) + (window.area - xi.area) * (n12 / (n12 + n22)))
}

#' @describeIn sae.v1ab.mean  The variance estimate
sae.v1ab.var <- function(xi.area, window.area, n11, n21, n12, n22){
  ahat11 <- (n11 / (n11 + n21))
  ahat12 <- (n12 / (n12 + n22))
  ssq11 <- (1 / (n11 + n21 -1)) * ( n11 * (1 - ahat11)^2 + n21 * (ahat11^2)) #sample variances
  ssq12 <- (1 / (n12 + n22 -1)) * ( n12 * (1 - ahat12)^2 + n22 * (ahat12^2)) #sample variances
  return( xi.area^2 * ssq11 + (window.area - xi.area)^2 * ssq12 )
}