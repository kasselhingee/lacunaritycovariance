#' @name placegrainsfromlib 

#' @title A function to help simulate Boolean models with user-provided grains
#' @aliases placegrainsfromlib
#' @export placegrainsfromlib
#' @author Kassel Liam Hingee

#' @description
#' A Boolean model has two components, a point process (called germs) and a process that creates
#'  independent identically distributed grains that are centred on the germs.
#' The point process of germs can be easily simulated using a variety of \code{spatstat} functions
#'  (note that this simulation window must include a buffer because grains centred outside an observation window can still be observed).
#'   The function here, \code{placegrainsfromlib},
#'   can then be used to randomly select grains from a library and place them around each point.
#'   The buffer zone must then be cropped out to get a simulation of the Boolean model in the desired observation window.


#' @param pp A point pattern (in \code{ppp} format).
#' @param grainlib A list of grains (in \code{\link[spatstat]{solist}} format) that grains will be selected from
#' @param replace passed directly to \code{\link[base]{sample}}. When TRUE grains are chosen from library with replacement.
#' @param prob A list of probability weights for each grain in the library. Passed directly to \code{\link[base]{sample}}.
#'  If NULL the grains are selected with equal probability.
#' @param w Optional desired observation window. If this is non-null then any grains with Frame outside the Frame of \code{w} will be ignored.
#' This reduces polygonal intersection calculations for very large buffer distances

#' @details \code{placegrainsfromlib} randomly samples from a library of grains (\code{grainlib}) and places these on the points in \code{pp}.

#' @return Returns an \code{owin} object.
#' @rdname placegrainsfromlib
#' 
#' @examples
#' #Generate a germ-grain models where germs are a Poisson point process
#' # and grains are 2 or 3 different disc sizes.
#' grainlib <- solist(disc(radius = 1), disc(radius = 1.9), disc(radius = 0.2))
#' bufferdist <- 2 #chosen to be larger than the largest radius in library
#' 
#' w <- owin(xrange = c(0, 10), yrange = c(0, 10))
#' 
#' #simulate the germ process in an enlarged window
#' pp <- rpoispp(lambda = 0.1, win = dilation(w, bufferdist), nsim = 1, drop = TRUE)
#'
#' plot(w)
#' xibuffer <- placegrainsfromlib(pp, grainlib)
#' plot(xibuffer, add = TRUE, lty = "dashed")
#' 
#' #get final simulation by intersection with desired window
#' xi <- intersect.owin(xibuffer, w)
#' plot(xi, hatch = TRUE, add = TRUE)
#' 
#' #demonstration that involves rasterisation.
#' xibuffer <- placegrainsfromlib(pp, grainlib, xy = as.mask(w, eps = 0.1))
#' plot(xibuffer)
#' plot(w, add = TRUE)
#' 
#' #Demo of covariance and set covariance computations: test on Boolean model
#' lambda <- 0.1
#' discr <- 10
#' weights <- c(0.9999, 0.0001)
#' grainlib <- solist(disc(r = discr), disc(r = 2*discr))
#' meangrainarea(grainlib, weights)
#' plot(meangrainsetcov(grainlib, weights, xy = as.mask(w, eps = 0.1)))
#' truecovartest <- grainlib.covar(lambda, grainlib, weights, xy = as.mask(w, eps = 0.1))
#' truecovariance <- bddcovar(
#'                    c(-10, 10), c(-10, 10), c(0.1, 0.1), lambda, discr)
#' plot(solist(truecovartest, truecovariance), clipwin = disc(r = 3))
#' plot(truecovartest - truecovariance, clipwin = disc(r = 3))
#' range(truecovartest - truecovariance)

#' @keywords spatial nonparametric datagen
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
meangrainarea <- function(grainlib, weights = rep(1/length(grainlib), length(grainlib))){
  grainareas <- vapply(grainlib, area.owin, 0.0)
  return(sum(grainareas * as.vector(weights)))
}

#' @describeIn placegrainsfromlib Compute the mean set covariance of the random grain given by the library.
#' xy is required because the set covariance function must rasterise the owin objects
meangrainsetcov <- function(grainlib, weights = rep(1/length(grainlib), length(grainlib)), xy){
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
grainlib.covar <- function(lambda, grainlib, weights, xy){
  p <- 1 - exp(- lambda * meangrainarea(grainlib, weights))
  meangsetcov <- meangrainsetcov(grainlib, weights, xy)
  racscov <- 2 * p - 1 + (1 - p)^2 * exp(lambda * meangsetcov) #this formula from 47 of Bohm 2004 Kernel Estimation of Stationary RACS paper
  return(racscov)
}

