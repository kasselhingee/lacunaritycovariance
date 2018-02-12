#' @title Variance Estimates for Observed Area
#' @export sae.all.mean sae.all.var
#' @description Estimates the variance of the area of a cover type observed in a thematic map created using a fallible classifier from remote sensing.
#' @author{Kassel Hingee}

#' @param xi An observation of the RACS of interest in 1, 0, or NA valued pixels.
#' Pixels must be square. Must be a spatstat im object.
#' @param obswin Observation window
#' @param corrrad Radius of the step function in the correlation (in the same units as xi)
#' @param corrstepheight Height of the step in the correlation
#' @param n11 The number of samples that were truly class 1 and fallibly classified as class 1
#' @param n21 The number of samples that were truly class 2 and fallibly classified as class 1
#' @param n12 The number of samples that were truly class 1 and fallibly classified as class 2
#' @param n22 The number of samples that were truly class 2 and fallibly classified as class 2
