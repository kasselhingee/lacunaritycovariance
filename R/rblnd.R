#' @title Simulate a Boolean model of discs with log normal disc radii
#' @export rblnd
#' @importFrom stats rlnorm

#' @description Simulates a Boolean model of discs with log normal radii by first simulating a Poisson point process and then placing discs
#'  of random radii around each point (the radii are generated using a log normal distribution).
#' @param  obswin An \code{owin} object specifying the desired simulation region
#' @param  bufferdist A distance to expand \code{obswin} so that discs with centres near \code{obswin} are also simulated.
#' @param lambda Intensity of the Poisson point process, passed to \link[spatstat.core]{rpoispp}.
#'   It could be either a single positive number, or any other object that \link[spatstat.core]{rpoispp} can understand.
#' @param  meanlog For the distribution of radii. The logarithm of the distribution is set to have mean \code{meanlog}.
#' @param  sdlog For the distribution of radii. The logarithm of the distribution is set to have standard deviation \code{sdlog}
#' @param seed Optional input (default is \code{NULL}). Is an integer passed to \code{\link[base]{set.seed}}. Used to reproduce patterns exactly.
#' 
#' @section Warning: A good choice of \code{bufferdist} is required and will be sensitive to the distribution of radii. 


#' @details The point process needs to be simulated in a larger region than the desired observation window to account for the possibility of discs that intersect the observation window, but have germs outside the observation window.
#'
#' The point process of germs is generated using spatstat's \code{\link[spatstat.core]{rpoispp}}.
#' @return Returns an \code{owin} object cropped to \code{obswin}.

#' @examples
#' w <- owin(xrange = c(0, 10), yrange = c(0, 10))
#' xi <- rblnd(w, 2, 0.3, -1, 0.2)


#' @keywords spatial datagen
rblnd <- function(obswin, bufferdist, lambda, meanlog, sdlog, seed = NULL){
  #have to simulate in a much larger area than the observation window (because grains with centres outside the window should still be observed)
  wsim <- Frame(dilation(obswin, bufferdist)) #i reckon faster to use rectangular region (the non-rectangular probably simulates in a rectangular region and then rejects anyway)
  if (!missing(seed)){set.seed(seed)}
  pp <- rpoispp(lambda, win = wsim, nsim = 1, drop = TRUE) #prepare a random radius for each point

  if (!missing(seed)){set.seed(seed)}
  radius <- rlnorm(pp$n, meanlog = meanlog, sdlog = sdlog) #prepare a random radius for each point

  #calculating grains
  pointlocations <- cbind(X = pp$x, Y = pp$y)
  pointlocations <- split(cbind(pointlocations), row(pointlocations)) #split matrix into a list of the rows
  grains <- mapply(disc, radius = radius, centre = pointlocations, SIMPLIFY = FALSE) #calculate grains with their locations

  #take union of all grains
  xisim <- union.owin(as.solist(grains))

  xi <- intersect.owin(xisim, obswin)
  return(xi)
}
