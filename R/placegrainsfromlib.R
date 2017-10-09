#' @name placegrainsfromlib 

#' @title A function to help simulate Boolean models with user-provided grains
#' @aliases placegrainsfromlib
#' @export placegrainsfromlib
#' @author Kassel Hingee

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
#' 

#' @keywords spatial nonparametric datagen
placegrainsfromlib <- function(pp, grainlib, replace = TRUE, prob = NULL){
  if (pp$n == 0){
    warning("there were no points in the point process - returning empty window")
    return(NULL)
  }
  grains <- sample(grainlib, size = pp$n, replace = replace, prob = prob)
  pointlocations <- cbind(X = pp$x, Y = pp$y)
  pointlocations <- split(cbind(pointlocations), row(pointlocations)) #split matrix into a list of the rows
  shiftedgrains <- as.solist(mapply(shift.owin, grains, vec = pointlocations, SIMPLIFY = FALSE))
  placedgrains <- union.owin(shiftedgrains)
# Note on union.owin: for pixel masks it uses inside.owin(xcol, yrow, A) | inside.owin(xcol,yrow,B) to determine union mask. It does this recursively.
# Inside owin uses a lot of checking about polygons etc.
  return(placedgrains)
}
