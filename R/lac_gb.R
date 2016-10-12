#' @title Gliding Box lacunarity from a black and white image
#' @export lacunarityGB
#'
#' @description Calculates the gliding box lacunarity
#'
#' @examples
data(balcattapark_coarse)
balcattapark_coarse$vegmask

Algorithm Plan:
+ do FFT (convolve.im) with a uniform kernel.
+ sum over the result in an eroded window
+ change kernel and repeat

+ get lidar data at some point


