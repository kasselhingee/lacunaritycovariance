#' @title Place grains randomly on a point pattern
#' @aliases placegrainsfromlib
#' @export placegrainsfromlib covar.grainlib meanarea.grainlib
#' @author Kassel Liam Hingee

#' @description
#' Places subsets (grains) of two dimension space randomly on a given point pattern.
#' This is useful for simulating germ-grain models such as Boolean models.
#' Also described here are functions for computing summary properties of the a list of grains.
#' 
#' 


#' @param pp A point pattern (in \code{ppp} format).
#' @param grainlib A list of grains as \code{owin} objects in a \code{\link[spatstat.geom]{solist}}. 
#' @param replace passed directly to \code{\link[base]{sample}}. When TRUE grains are chosen from library with replacement.
#' @param prob A list of probability weights for each grain in the library. Passed directly to \code{\link[base]{sample}}.
#'  If NULL the grains are selected with equal probability.
#' @param w Optional desired observation window. If this is non-null then any grains with Frame outside the Frame of \code{w} will be ignored.
#' This reduces polygonal intersection calculations for very large buffer distances
#' @param xy An \code{im} or binary mask object that is used to specify the pixel array of objects.
#' @param weights Probability of selecting each grain in the library
#' @param lambda Intensity of germs of a Boolean model - for computing the covariance of a Boolean model that has grain distribution given by \code{grainlib} and \code{weights}.

#' @details 
#' Germ-grain models have two components, a point process (called germs) and a process that creates
#'  grains that are centred on the germs.
#' The point process of germs can be easily simulated using a number of \pkg{spatstat} functions 
#' (e.g. \code{\link[spatstat.core]{rpoispp}} for Boolean models).
#' To simulate a germ-grain model in a window \eqn{W} the germ process must be simulated in a larger window 
#' because grains centred outside \eqn{W} can intersect \eqn{W}.
#' The result must then be cropped to \eqn{W} to achieve a realisation of the germ-grain process within \eqn{W}.
#' 
#' \code{placegrainsfromlib} randomly samples from a library of grains (\code{grainlib}) and places these on the points in \code{pp}.
#' Sampling of the grain is independent of the location of the point in \code{pp}.
#' It can be used to simulate the grain process in some germ-grain models.
#' 
#' 
#' @return Returns an \code{owin} object.
#' 
#' @examples
#' # Simulate a germ-grain model where germs are a Poisson point process
#' # and grains are randomly selected from 3 different disc sizes.
#' grainlib <- solist(disc(radius = 1), disc(radius = 1.9), disc(radius = 0.2))
#' bufferdist <- 2 #chosen to be larger than the largest radius in library
#' w <- owin(xrange = c(0, 10), yrange = c(0, 10))
#' 
#' # Simulate the germ process in the window plus a buffer region around window
#' pp <- rpoispp(lambda = 0.1, win = dilation(w, bufferdist), nsim = 1, drop = TRUE)
#' xi_withbuffer <- placegrainsfromlib(pp, grainlib)
#' 
#' # Simulation of germ-grain model is the part within the window
#' xi <- intersect.owin(xi_withbuffer, w)
#'
#' # Computation of properties from parameters 
#' lambda <- 0.1
#' discr <- 10
#' weights <- c(0.9999, 0.0001)
#' grainlib <- solist(disc(r = discr), disc(r = 2*discr))
#' meanarea.grainlib(grainlib, weights)
#' truecovar <- covar.grainlib(lambda, grainlib, weights, xy = as.mask(w, eps = 2))

#' @keywords spatial nonparametric datagen
#' @describeIn placegrainsfromlib Place grains randomly from a list of grains.
placegrainsfromlib <- function(pp, grainlib,
                               replace = TRUE, prob = NULL, w = NULL, xy = NULL){
  if (pp$n == 0){
    warning("there were no points in the point process - returning empty window")
    return(NULL)
  }
  if (!is.null(xy)){
    grainlib <- solapply(grainlib, tocompatiblepixelarray, xy = xy)
  }
  grains <- sample(grainlib, size = pp$n, replace = replace, prob = prob)
  pointlocations <- coords(pp)
  if (!is.null(xy)){
    pointlocations$x <- vapply(pointlocations$x, function(inx) xy$xcol[which.min(abs(inx - xy$xcol))], 0.00)
    pointlocations$y <- vapply(pointlocations$y, function(inx) xy$yrow[which.min(abs(inx - xy$yrow))], 0.00)
  }
  pointlocations <- split(pointlocations, 1:nrow(pointlocations)) #split matrix into a list of the rows
  shiftedgrains <- as.solist(mapply(shift.owin, grains, vec = pointlocations, SIMPLIFY = FALSE))
  if (!is.null(w)){shiftedgrains <- shiftedgrains[vapply(shiftedgrains, isinwindowbbox, FUN.VALUE = c(FALSE), w = w)]} #this line removes grains that aren't likely to intersect window
  placedgrains <- union.owin(shiftedgrains)
# Note on union.owin: for pixel masks it uses inside.owin(xcol, yrow, A) | inside.owin(xcol,yrow,B) to determine union mask. It does this recursively.
# Inside owin uses a lot of checking about polygons etc.
  return(placedgrains)
}

#tests if in window bounding box using rectangle arithmetic to avoid excessive computation
isinwindowbbox <- function(grain, w){
  intersection <- intersect.owin(Frame(w), Frame(grain), fatal = FALSE)
  return(!is.null(intersection))
}


#convert grain to raster on the same sized pixel grid and origin
tocompatiblepixelarray <- function(grain, xy){
  stopifnot(is.owin(grain))
  grainw <- Frame(grain)
  xcol <- seq(floor(grainw$xrange[[1]]/xy$xstep)*xy$xstep,
              ceiling(grainw$xrange[[2]]/xy$xstep)*xy$xstep, by = xy$xstep)
  yrow <- seq(floor(grainw$yrange[[1]]/xy$ystep)*xy$ystep,
              ceiling(grainw$yrange[[2]]/xy$ystep)*xy$ystep, by = xy$ystep)
  return(as.mask(grain, xy = list(x = xcol, y = yrow)))
}

#' @describeIn placegrainsfromlib Compute mean area of a random grain given by the library
meanarea.grainlib <- function(grainlib, weights = rep(1/length(grainlib), length(grainlib))){
  grainareas <- vapply(grainlib, area.owin, 0.0)
  return(sum(grainareas * as.vector(weights)))
}

#' @describeIn placegrainsfromlib Computes the mean of the set covariance of the grains in \code{grainlib}.
#' \code{xy} is required because the set covariance function must rasterise the \code{owin} objects.
meansetcov.grainlib <- function(grainlib, weights = rep(1/length(grainlib), length(grainlib)), xy){
  grainlib <- solapply(grainlib, tocompatiblepixelarray, xy = xy)
  setcovagrain <- solapply(grainlib, setcov)
  setcovagrain <- as.solist(mapply(function(x, y) x * y, setcovagrain, weights, SIMPLIFY = FALSE))
  meansetcov <- 0 * as.im(tocompatiblepixelarray(do.call(union.owin,lapply(setcovagrain, Frame)), xy = xy)) #blank mean setcov
  setcovagrain <- solapply(setcovagrain, as.im, xy = meansetcov, na.replace = 0) #pad with zeros
  for (i in 1:length(grainlib)){
    meansetcov <- meansetcov + setcovagrain[[i]]
  }
  return(meansetcov)
}


#' @describeIn placegrainsfromlib Compute the covariance of a Boolean model with random grain given by the library
covar.grainlib <- function(lambda, grainlib, weights, xy){
  p <- 1 - exp(- lambda * meanarea.grainlib(grainlib, weights))
  meangsetcov <- meansetcov.grainlib(grainlib, weights, xy)
  racscov <- 2 * p - 1 + (1 - p)^2 * exp(lambda * meangsetcov) #this formula from 47 of Bohm 2004 Kernel Estimation of Stationary RACS paper. Holds for any stationary Boolean model and any grain distribution - including ones that are anisotropic
  return(racscov)
}

