#' @title Simulate Boolean Model with Grains Scaled According to a Truncated
#'   Pareto Distribution
#' @export rbpto bpto.coverageprob bpto.covar bpto.germintensity
#' @description Functions for simulation and computing theoretical values of a
#'   Boolean model with identically shaped grains with size given by a
#'   truncated Pareto distribution.
#'
#' @param lambda Intensity of the germ process (which is a Poisson point
#'   process)
#' @param grain A single \code{owin} object that gives the shape and size of the grain
#'   at scale 1
#' @param xm A parameter governing the shape of the Pareto distribution used -
#'   see details
#' @param alpha A parameter governing the shape of the Pareto distribution used
#'   - see details
#' @param lengthscales A list of scales of the \code{grain} for which to
#'   approximate the Pareto distribution: The grain for a germ is chosen by
#'   selecting a scaled version of \code{grain} where \code{lengthscales}
#'   specifies the possible scales and the Pareto distribution is used to
#'   specify the probability of selection of each scale.
#' @param coverp Coverage probability of the Boolean model.
#' @param win The window to simulate in (an \code{owin} object)
#' @param seed Optional input (default in NULL). Is an integer passed to
#'   \code{\link[base]{set.seed}}. Used to reproduce patterns exactly.
#' @param xy A raster object that specifies pixel coordinates of the final
#'   simulated binary map. It is used the same way as \code{xy} is
#'   \code{\link[spatstat.geom]{as.mask}} in \pkg{spatstat}. If non-null then the
#'   computations will be performed using rasters. Otherwise if \code{grain} and
#'   \code{win} are polygonal then computations may be all polygonal.

#' @details
#' The parameters \code{xm} and \code{alpha} are such that the CDF of the Pareto distribution is \eqn{P(s <= x) = 1 - (xm / x)^{alpha}}.
#' The distribution of grains scales is a step-function approximation to the CDF with steps at \code{lengthscales}.
#' 
#'
#' 
#' @return 
#' An \code{owin} object.
#' @keywords spatial datagen


#' @examples
#' lambda <- 0.2
#' win <- square(r = 10)
#' grain <- disc(r = 0.2)
#' xm <- 0.01
#' alpha <- 2
#' lengthscales <- seq(1, 5, by = 0.1)
#' xi <- rbpto(lambda, grain, win, xm, alpha, lengthscales = lengthscales)
#' 
#' # Compute properties of the Boolean model from parameters
#' bpto.coverageprob(lambda, grain, xm, alpha, lengthscales = lengthscales)
#' covar <- bpto.covar(lambda, grain, xm, alpha, lengthscales = lengthscales,
#'                     xy = as.mask(win, eps = 2))


#' @describeIn rbpto Simulate Boolean model with grain size distributed according to a truncated Pareto distribution.
rbpto <- function(lambda, grain, win, xm, alpha, lengthscales,
                  seed = NULL, xy = NULL){
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

#' @describeIn rbpto  The coverage probability of the Boolean model with grain size distributed according to a truncated Pareto distribution.
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

#' @describeIn rbpto  The germ intensity of the Boolean model with grain size distributed according to a truncated Pareto distribution.
bpto.germintensity <- function(coverp, grain, xm, alpha,
                              lengthscales = 1:500){
  #first get mean area. Need weight for each discrete grain.
  #check that smallest scale is larger than xm
  stopifnot(lengthscales[1] >= xm)
  
  #get weights of grain sizes from pmf
  weights <- alpha * xm ^ alpha / (lengthscales ^ (alpha + 1) )
  weights <- weights / sum(weights) #standardise
  
  meangrainarea <- sum(lengthscales * lengthscales * area.owin(grain) * weights)
  #coverp formula: p = 1 - exp(- lambda * meangrainarea)
  lambda <-   - 1 * log(1 - coverp)/ meangrainarea
  return(lambda)
}

#' @describeIn rbpto  The covariance of the Boolean model with grain size distributed according to a truncated Pareto distribution.
#' \code{xy} is required to specify resolution and offset of pixel grid.
bpto.covar <- function(lambda, grain, xm, alpha, lengthscales = 1:500, xy){
  #check that smallest scale is larger than xm
  stopifnot(lengthscales[1] >= xm)
  
  #first get grainlib with weights
  #get scaled versions of grain
  grainlib <- mapply(scalardilate, X = list(grain), f = lengthscales, SIMPLIFY = FALSE)
  grainlib <- as.solist(grainlib)
  
  #get weights of these grains from pmf
  weights <- alpha * xm ^ alpha / (lengthscales ^ (alpha + 1) )
  weights <- weights / sum(weights) #standardise
  
  covar <- covar.grainlib(lambda, grainlib, weights, xy)
  return(covar)
}
