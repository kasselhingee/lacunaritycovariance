#' @name rgermgrain
NULL

#' @title Functions for flexibly generating a simulation of germ grain model
#' @aliases placegrainsfromlib
#' @export placegrainsfromlib
#' @author Kassel Hingee

#' @description
#' Use point process simulations from spatstat to generate a point process in desired region (must include a buffer). \code{placegrainsfromlib} can then be used to randomly select grains from a library and place them around each point. Crop out the buffer zone to get a simulation.


#' @param pp A point pattern (in \code{ppp} format).
#' @param grainlib A list (aka. library/population) of grains (in \code{solist} format) that grains will be selected from
#' @param replace passed directly to \code{\link[base]{sample}}. When TRUE grains are chosen from library with replacement.
#' @param prob A list of probability weights for each grain in the library. Passed directly to \code{\link[base]{sample}}. If NULL I'm pretty sure grains all have equal probability.

#' @details \code{placegrainsfromlib} randomly samples from a library of grains (\code{grainlib}) and places these on the points in \code{pp}.

#' @return Returns an \code{owin} object.
#' @rdname rgermgrain
#' 
#' @examples
#' #Generating a germ-grain models where germs are a Poisson Point process, and grains are 2 or 3 different disc sizes.
#' grainlib <- solist(disc(radius=1),disc(radius=1.9),disc(radius=0.2))
#' bufferdist <- 2 #chosen to be larger than the largest radius in library
#' 
#' w <- owin(xrange=c(0,10),yrange=c(0,10))
#' 
#' pp <- rpoispp(lambda=0.1,win=dilation(w,bufferdist),nsim=1,drop=TRUE)
#' xibuffer <- placegrainsfromlib(pp,grainlib)
#' xi <- intersect.owin(xibuffer,w)
#' 
#' plot(w)
#' plot(xibuffer,add=TRUE)
#' plot(xi,hatch=TRUE,add=TRUE)
#' plot(pp,pch="+",add=TRUE)

#' @keywords spatial nonparametric 
placegrainsfromlib <- function(pp,grainlib,replace=TRUE,prob=NULL){
  if (pp$n == 0){
    warning("there were no points in the point process - returning empty window")
    return(NULL)
  }
  grains <- sample(grainlib,size=pp$n,replace=replace,prob=prob) 
  pointlocations <- cbind(X=pp$x,Y=pp$y)
  pointlocations <- split(cbind(pointlocations),row(pointlocations)) #split matrix into a list of the rows
  shiftedgrains <- as.solist(mapply(shift.owin,grains,vec=pointlocations,SIMPLIFY=FALSE))
  placedgrains<- union.owin(shiftedgrains)  
  return(placedgrains)
}

# Note on union.owin: for pixel masks it uses inside.owin(xcol, yrow, A) | inside.owin(xcol,yrow,B) to determine union mask. It does this recursively.
# Inside owin uses a lot of checking about polygons etc.


