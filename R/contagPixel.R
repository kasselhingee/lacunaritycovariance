#' @title Pixel Contagion
#' @export contagpixelgrid adjacency
#' 
#' @description Function for calculating the standard LPI contagion from a pixel grid.
#' Refer to FRAGSTATS for a full description of contagion.
#' 
#' @param xi an owin mask representing a realisation of RACS
#' @param w the window of observation
#' @param normalise If \code{TRUE} will divide result by 2*log(2) and add 1 to make contagion between 0 and 1 for any binary map

#' @details The `double count' method as described in the FRAGSTATS manual on adjacency matrix.
#' Considering not xi as another phase and using a 4-neighbourhood.
#' @examples
#' xi <- heather$coarse
#' w <- owin(xrange = c(0,7),yrange=c(0,16))
#' adjmat <- adjacency(xi,w)
#' adjmat
#' contagion <- contagpixelgrid(xi,w)
#' contagion

#' #Finer resolution
#' xi <- heather$medium
#' w <- owin(xrange = c(0,7),yrange=c(0,16))
#' contagion <- contagpixelgrid(xi,w)
#' contagion

#' @section Warning: Will fail if there are no adjacencies
contagpixelgrid <- function(xi, w, normalise=FALSE){
  stopifnot(is.mask(xi))
  out <- harmonise(xi,w)
  xi <- out[[1]]
  w <- out[[2]]
  adjmat <- adjacency(xi,w)
  # num pixels in w?
  propOfXi <- sum(as.matrix(intersect.owin(xi,w)))/sum(as.matrix(w))
  propOfNotXi <- sum(as.matrix(intersect.owin(complement.owin(xi),w)))/sum(as.matrix(w))
  
  contag <- 0
  #xi with xi part of contagion
  contag <- contag + sum(propOfXi*adjmat[1,]/rowSums(adjmat)[1]*log(propOfXi*adjmat[1,]/rowSums(adjmat)[1]))
  contag <- contag + sum(propOfNotXi*adjmat[2,]/rowSums(adjmat)[2]*log(propOfNotXi*adjmat[2,]/rowSums(adjmat)[2]))
  if (normalise) {contag <- 1+contag/(2*log(2))}
  return(contag)
}


#' @describeIn contagpixelgrid Calculates the adjacency matrix used in the pixel contagion
adjacency <- function(xi,w){
  stopifnot(is.mask(xi))
  xic <- intersect.owin(complement.owin(xi),w)
  xi <- intersect.owin(xi,w)
  #neighbours of points in xi that are also in xi, in each direction
  #(4 neighbourhood)
  numNnbr <- sum(as.matrix(intersect.owin(xi,shift.owin(xi,vec=c(0,-xi$ystep)))))
  numSnbr <- sum(as.matrix(intersect.owin(xi,shift.owin(xi,vec=c(0,xi$ystep)))))
  numWnbr <- sum(as.matrix(intersect.owin(xi,shift.owin(xi,vec=c(xi$xstep,0)))))
  numEnbr <- sum(as.matrix(intersect.owin(xi,shift.owin(xi,vec=c(-xi$xstep,0)))))
  xinbrsofxi <- sum(c(numNnbr,numSnbr,numWnbr,numEnbr))  
  
  #neighbours of points in xi that are not in xi
  numNnbr <- sum(as.matrix(intersect.owin(xi,shift.owin(xic,vec=c(0,-xi$ystep)))))
  numSnbr <- sum(as.matrix(intersect.owin(xi,shift.owin(xic,vec=c(0,xi$ystep)))))
  numWnbr <- sum(as.matrix(intersect.owin(xi,shift.owin(xic,vec=c(xi$xstep,0)))))
  numEnbr <- sum(as.matrix(intersect.owin(xi,shift.owin(xic,vec=c(-xi$xstep,0)))))
  notxinbrsofxi <- sum(c(numNnbr,numSnbr,numWnbr,numEnbr)) 
  
  #nbrs of not xi points that are IN xi
  numNnbr <- sum(as.matrix(intersect.owin(xic,shift.owin(xi,vec=c(0,-xi$ystep)))))
  numSnbr <- sum(as.matrix(intersect.owin(xic,shift.owin(xi,vec=c(0,xi$ystep)))))
  numWnbr <- sum(as.matrix(intersect.owin(xic,shift.owin(xi,vec=c(xi$xstep,0)))))
  numEnbr <- sum(as.matrix(intersect.owin(xic,shift.owin(xi,vec=c(-xi$xstep,0)))))
  xinbrsofnotxi <- sum(c(numNnbr,numSnbr,numWnbr,numEnbr)) 
  
  #nbrs of not xi points that are NOT IN xi
  numNnbr <- sum(as.matrix(intersect.owin(xic,shift.owin(xic,vec=c(0,-xi$ystep)))))
  numSnbr <- sum(as.matrix(intersect.owin(xic,shift.owin(xic,vec=c(0,xi$ystep)))))
  numWnbr <- sum(as.matrix(intersect.owin(xic,shift.owin(xic,vec=c(xi$xstep,0)))))
  numEnbr <- sum(as.matrix(intersect.owin(xic,shift.owin(xic,vec=c(-xi$xstep,0)))))
  notxinbrsofnotxi <- sum(c(numNnbr,numSnbr,numWnbr,numEnbr)) 
  
  
  adjacencymat <- matrix(c(xinbrsofxi,notxinbrsofxi,xinbrsofnotxi,notxinbrsofnotxi), nrow=2, ncol=2, byrow = TRUE)
  colnames(adjacencymat) <- c("Xi","Not Xi")
  rownames(adjacencymat) <- c("Xi", "Not Xi")
  
  return(adjacencymat)
}

