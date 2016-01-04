#' @title Simulate a Boolean model of discs with log normal disc radii and Poisson germs
#' @description Simulates a Poisson point process and places discs of random radii around each point.  The radii are generated using a log normal distribution.
#' @param  window An \code{owin} object specifying the desired simulation region
#' @param  bufferdist A distance to expand \code{window} so that discs with centres near \code{window} are also simulated.
#' @param lambda Intensity of the Poisson point process, passed to \code{\link{spatstat}{rpoispp}}.  It could be either a single positive number, a function (x,y,...) or a pixel image (they must be defined for the dilated window)
#' @param  meanlog For the distribution of radii. The logarithm of the distribution has mean equal to \code{meanlog}.
#' @param  sdlog For the distribution of radii. The logarithm of the distribution has standard deviation equal to \code{sdlog}
#' @section Warning: A good choice of bufferdist is required and it probably depends on the distribution of radii. 


#' @details The point process needs to be simulated in a larger region than the desired window to account for the possibility of discs that intersect the window, but have germs outside the window.
#'
#' The point process is generated using spatstat's \code{\link[spatstat]{rpoispp}}.
#' @return Returns an \code{owin} object cropped to \code{window}.

#' @export rboollognormdiscs
rboollognormdiscs <- function(window,bufferdist,lambda,meanlog,sdlog){
  #have to simulate in a much larger area than the observation window (because grains with centres outside the window should still be observed)
  wsim <- Frame(dilation(window,bufferdist)) #i reckon faster to use rectangular region (the non-rectangular probably simulates in a rectangular region and then rejects anyway)
  pp <- rpoispp(lambda,win=wsim,nsim=1,drop=TRUE) #prepare a random radius for each point
 
  radius <- rlnorm(pp$n,meanlog=meanlog,sdlog=sdlog) #prepare a random radius for each point
   
  #calculating grains
  pointlocations <- cbind(X=pp$x,Y=pp$y)
  pointlocations <- split(cbind(pointlocations),row(pointlocations)) #split matrix into a list of the rows
  grains <- mapply(disc,radius = radius,centre=pointlocations,SIMPLIFY=FALSE) #calculate grains with their locations
  
  #take union of all grains
  xisim <- union.owin(as.solist(grains))
  
  xi <- intersect.owin(xisim,window)
  return(xi) 
}


#' @references [1] Molchanov, I. (1997) Statistics of the Boolean Model for Practitioners and Mathematicians. Wiley. (Chapter 8)
#' 
#' [2] **seminal paper - haven't read** Bindrich U. and Stoyan D. (1991) Stereology for pores in white bread: statistical analyses for the Boolean model by serial sections. J. Microscopy 162:231-239

#' @examples
#' w <- owin(xrange=c(0,10),yrange=c(0,10))
#' xi <- rboollognormdiscs(w,2,1,-1,0.5)
#' 
#' plot(w)
#' plot(xi,add=TRUE)


#' @keywords  spatial nonparametric 
