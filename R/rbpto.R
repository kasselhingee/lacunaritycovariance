#' @title Simulation of Boolean Model of Grains Scaled According to a Pareto Distribution
#' @export rbpto bpto.coverageprob
#' 
#' @param lambda Intensity of the germ process (which is a Poisson point process)
#' @param grain A single owin object that gives the shape and size of the grain at scale 1
#' @param win The window to simulate in (an owin object)
#' @param seed Optional input (default in NULL). Is an integer passed to \code{\link{base}{set.seed}}. Used to reproduce patterns exactly.
#' @param xy A raster object that specifies pixel coordinates of the final simulated binary map.
#' It is used the same way as \code{xy} is \code{\link{spatstat}{as.mask}} in spatstat.
#' If non-null then the computations will be performed using rasters. Otherwise if \code{grain} and \code{win} are polygonal then computations may be all polygonal.

#' @details
#' The parameters xm and alpha are such that the CDF of the Pareto distribuion is \eqn{P(s <= x) = 1 - (xm / x)^{alpha}}.
#' Approximates the grain distribution using the scales given by \code{lengthscales} and weighted by the probability distribution function of the Pareto distribution.

#' @return 
#' An owin object.


#' @examples
#' lambda <- 1
#' win <- square(r = 10)
#' #grain <- owin(xrange = c(-0.2, 0.2), yrange = c(-0.2, 0.2))
#' grain <- disc(r = 0.2)
#' xm <- 0.01
#' alpha <- 2
#' 
#' #system.time(xi <- rbpto(lambda, grain, win, xm, alpha, lengthscales = 1:100, xy = as.mask(win, eps = 0.1)))
#' system.time(xi <- rbpto(lambda, grain, win, xm, alpha, lengthscales = 1:100))
#' plot(xi)
#' 
#' bpto.coverageprob(lambda, grain, xm, alpha, lengthscales = 1:100)

rbpto <- function(lambda, grain, win, xm, alpha,
                  seed = NULL, xy = NULL, lengthscales = 1:500){
  #check that smallest scale is larger than xm
  stopifnot(lengthscales[1] >= xm)
  
  #get scaled versions of grain
  grainlib <- mapply(scalardilate, X = list(grain), f = lengthscales, SIMPLIFY = FALSE)
  grainlib <- as.solist(grainlib)
  
  #get weights of these grains from pmf
  weights <- alpha * xm ^ alpha / (lengthscales ^ (alpha + 1) )
  weights <- weights / sum(weights) #standardise
  
  #now simulate Boolean model!
  bufferdist <- diameter.owin(grain) * max(lengthscales)
  set.seed(seed)
  pp <- rpoispp(lambda, win = dilation.owin(win, bufferdist), nsim = 1, drop = TRUE)
  #plot(pp)
  #plot(add = TRUE, w)
  
  #place grains
  set.seed(seed) #this is not best way, must be a way to continue using sequence of pseudo-independent random numbers in already selected seed.
  xibuffer  <- placegrainsfromlib(pp, grainlib, prob = weights, w = win, xy = xy)
  #plot(w)
  #plot(xibuffer, add = TRUE, col = "black")
  #plot(xibuffer[square(r = 50)], col = "black")
  xi <- intersect.owin(xibuffer, win)
  return(xi)
}

#' #describeIn rbpto  The coverage probability of the Boolean model with scaled grains distributed according to Pareto distribution. Uses approximation of truncated length scales.
bpto.coverageprob <- function(lambda, grain, xm, alpha,
                              lengthscales = 1:500){
  #first get mean area. Need weight for each discrete grain.
  #check that smallest scale is larger than xm
  stopifnot(lengthscales[1] >= xm)
  
  #get weights of grain sizes from pmf
  weights <- alpha * xm ^ alpha / (lengthscales ^ (alpha + 1) )
  weights <- weights / sum(weights) #standardise
  
  meangrainarea <- sum(lengthscales * lengthscales * area.owin(grain) * weights)
  return(1 - exp(- lambda * meangrainarea))
}