#' @title Two Point Probability Contagion
#' @export contagTwoPtProb
#' @description Calculates the two point probability version of Contagion. 
#' @details Contagion is **. 
#' In this case the Pij are the probability that two points separated by a vector \code{v}
#' are in class \eqn{i} and class \eqn{j} respectively.
#' For an observation of a single phase \eqn{\Xi} this is easy to calculate.
#'  Class 1 corresponds to being in \eqn{\Xi}, class 0 to being outside \eqn{\Xi} 
#'  
#'  If \code{p} is unspecified then the estimate covariance at the origin, \eqn{C(o)}, is used.
#'  
#'  If \code{v} is unspecified then a map of contagion for many possible vectors is provided.
#'  
#'  If \code{normalise} is \code{TRUE} then result is divided by \eqn{-log(1/4)} and translated by 1 to force contagion 
#'  between 0 and 1.
#'  
#' @section Warning: there might still be some instability for covariance very close, but less than p
#
#' @param covariance is a map of covariance in spatstat \code{im} format
#' @param p is an estimated coverage fraction. If none is provided an estimate is made using covariance (see details).
#' @param v is an optional input. It is a vector in format c(x,y)
#' @param normalise If true normalises the results so that all RACS return a value between 0 and 1
#' @return If \code{v} is included then returns a single number, otherwise a map of the two point contagion
#' 
#' @seealso \code{\link{contagdiscstate}} 

#' @examples 
#' xi <- heather$coarse
#' covariance <- racscovariance(xi,Frame(xi))
#' twoptcontagion <- contagTwoPtProb(covariance)
#' p <- coveragefrac(xi,Frame(xi))
#' twoptcontagion <- contagTwoPtProb(covariance,p)
#' # plot(twoptcontagion)
#' # plot(twoptcontagion,clipwin=owin(xrange=c(-0.5,0.5),yrange=c(-0.5,0.5)),main="zoom")
#' v <- c(5,15)
#' twoptcontagion[ppp(v[1],v[2],window=owin(xrange=c(v[1]-1,v[1]+1),yrange=c(v[2]-1,v[2]+1)))]
#' 
#' contagTwoPtProb(covariance,p=p, v=c(5,15))
#' #result for both should be -0.9985666
#' 
contagTwoPtProb <- function(covariance,p=NULL,v=NULL,normalise=FALSE){
  if (is.null(p)){p <- covariance[ppp(0,0)]}
  if (is.null(v)){
    originOnlyInXi <- covariance
    originOnlyInXi[p-covariance>0] <- (p-covariance[p-covariance>0])*log(p-covariance[p-covariance>0])
    originOnlyInXi[p-covariance <= 0] <- 0
    #vOnlyInXi equals originOnlyInXi
    neitherInXi <- 1-2*p + covariance
    neitherInXi <- covariance
    neitherInXi[1-2*p + covariance>0] <- (1-2*p + covariance[1-2*p + covariance>0])*log(1-2*p + covariance[1-2*p + covariance>0])
    neitherInXi[1-2*p + covariance <= 0] <- 0
    
    unnormalisedContag <- eval.im(covariance*log(covariance)+
      2*originOnlyInXi +
        neitherInXi)
    
    if (normalise){return(1+1/(log(1/4)) *unnormalisedContag)}
    return(unnormalisedContag)
  }
  else {
    bothInXiProb <- covariance[
      ppp(v[1],v[2],window=owin(xrange=c(v[1]-1,v[1]+1),yrange=c(v[2]-1,v[2]+1)))]
    if (p < bothInXiProb){#this is hopefully just a machine calculation error
      warning("coverage fraction estimate is smaller than the covariance")
      return(p*log(p)+(1-p)*log(1-p))
    }
    originOnlyInXi <- p-bothInXiProb
    vOnlyInXi <- p-bothInXiProb
    neitherInXi <- 1-2*p + bothInXiProb
    
    unnormalisedContag <- bothInXiProb*log(bothInXiProb)+
      originOnlyInXi*log(originOnlyInXi)+ 
      vOnlyInXi*log(vOnlyInXi)+
      neitherInXi*log(neitherInXi)
    if (normalise) {return(1+1/(log(1/4)) *unnormalisedContag)}
    else{return(unnormalisedContag)}
  }    
}

# my arithmetic tells me it should simplify to
# 1/(-log(1/4)) * (
#   (1+bothInXiProb-2*p)*log(1+bothInXiProb-2*p)
#   + 2*(p-bothInXiProb)*log(p-bothInXiProb)
#   + bothInXiProb*log(bothInXiProb)
# )
#and it matches :)


