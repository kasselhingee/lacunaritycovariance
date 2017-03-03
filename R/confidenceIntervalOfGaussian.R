#' @title Confidence intervals for Guassian Distributions 
 
#' @description A simple function for getting symmetric confidence intervals for a given confidence level and standard deviation of a Gaussian distribution

#' @param confidenceLevel The desired confidence level (between 0 and 1)
#' @param sd The standard deviation of the Guassian distribution.
#' @return A list of length two: (lower bound, upper bound)
#' @examples
#' confidenceIntervalOfGaussian(0.95,1)
#' 
#' @export confidenceIntervalOfGaussian
#' @importFrom stats qnorm
 confidenceIntervalOfGaussian <- function(confidenceLevel,sd){
   minBound <- qnorm((1-confidenceLevel)/2,sd=sd)
   maxBound <- -minBound
   return(c(minBound,maxBound))
 }
