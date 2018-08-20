#' @title Pixel Adjacency Contagion
#' @export contagpixelgrid adjacency
#' 
#' @description Function for calculating the traditional contagion LPI from a binary map [cite fragstats manual].
#' 
#' @param xi A binary map of an observation of a RACS of interest. See
#'   \code{\link{stationaryracsinference-package}} for details.
#'   If \code{xi} is an owin object it must be of \code{mask} type.
#' @param obswin If \code{xi} is an \code{owin} object then \code{obswin} is an
#'   \code{owin} object that specifies the observation window.
#' @param normalise If \code{TRUE} will divide result by \eqn{2 ln(2)} and add 1 to make contagion between 0 and 1 for any binary map

#' @details The unnormalised contagion LPI of categorical map is defined as
#' \deqn{\sum_i \sum_j Pij ln(Pij),} where \eqn{Pij} is the probability of
#' randomly selected adjacent pixels being in class \eqn{i} and class \eqn{j}
#' respectively, and \eqn{m} is the number of classes.
#'
#' Here \eqn{m = 2} as \eqn{xi} is a binary map and we have defined 'adjacent'
#' pixels using the 4-neighbourhood regime.
#' 
#' Contagion is calculated from an adjacency matrix created using \code{adjacency}.
#' 
#' **What is the `double count' method as described in the
#' FRAGSTATS manual on adjacency matrix? 

#' @examples
#' xi <- heather$coarse
#' obswin <- owin(xrange = c(0,7),yrange=c(0,16))
#' adjmat <- adjacency(xi,obswin)
#' adjmat
#' contagion <- contagpixelgrid(xi,obswin)
#' contagion

#' #Finer resolution
#' xi <- heather$medium
#' obswin <- owin(xrange = c(0,7),yrange=c(0,16))
#' contagion <- contagpixelgrid(xi,obswin)
#' contagion

#' @section Warning: Will fail if there are no adjacencies
#' @describeIn contagpixelgrid Traditional contagion LPI of a binary map.
contagpixelgrid <- function(xi, obswin, normalise=FALSE){
  stopifnot(is.mask(xi))
  out <- harmonise(xi,obswin)
  xi <- out[[1]]
  obswin <- out[[2]]
  adjmat <- adjacency(xi,obswin)
  # num pixels in obswin?
  propOfXi <- sum(as.matrix(intersect.owin(xi,obswin)))/sum(as.matrix(obswin))
  propOfNotXi <- sum(as.matrix(intersect.owin(complement.owin(xi),obswin)))/sum(as.matrix(obswin))
  
  contag <- 0
  #xi with xi part of contagion
  contag <- contag + sum(propOfXi*adjmat[1,]/rowSums(adjmat)[1]*log(propOfXi*adjmat[1,]/rowSums(adjmat)[1]))
  contag <- contag + sum(propOfNotXi*adjmat[2,]/rowSums(adjmat)[2]*log(propOfNotXi*adjmat[2,]/rowSums(adjmat)[2]))
  if (normalise) {contag <- 1+contag/(2*log(2))}
  return(contag)
}


#' @describeIn contagpixelgrid Calculates the adjacency matrix used in the pixel contagion
adjacency <- function(xi, obswin = NULL){
  if("im" %in% class(xi)){isbinarymap(xi, requiretrue = TRUE)}
  else if (is.owin(xi) && is.null(obswin)){stop("obswin must be included if xi is an owin object.")}
  
  if (is.owin(xi){
    xi <- as.im(xi, value = TRUE, na.replace = FALSE)
    xi[setminus.owin(Frame(xi), obswin)] <- NA
  }
  xic <- intersect.owin(complement.owin(xi),obswin)
  xi <- intersect.owin(xi,obswin)
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

