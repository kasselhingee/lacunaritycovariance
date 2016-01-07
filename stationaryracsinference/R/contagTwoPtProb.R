#' @title Two Point Probability Contagion
#' @export contagTwoPrProb
#' @description Calculated the two point probability version of Contagion at a single vector. 
#' @details Contagion is **. 
#' In this case the Pij are the probability that two points separated by a vector \code{v}
#' are in class \eqn{i} and class \eqn{j} respectively.
#' For an observation of a single phase \eqn{\Xi} this is easy to calculate.
#'  Class 1 corresponds to being in \eqn{\Xi}, class 0 to being outside \eqn{\Xi} 
#' 
#
#' @param v is a vector in format c(x,y)
#' @param covariance is a map of covariance in spatstat \code{im} format
#' @param p is an estimated coverage fraction. An estimate could also use \eqn{C(o)}.
#' @return a single number
contagTwoPtProb <- function(v,covariance,p){
  bothInXiProb <- covariance[
    ppp(v[1],v[2],window=owin(xrange=c(v[1]-1,v[1]+1),yrange=c(v[2]-1,v[2]+1)))]
  originOnlyInXi <- p-bothInXiProb
  vOnlyInXi <- p-bothInXiProb
  neitherInXi <- 1-2*p + bothInXiProb
  
  unnormalisedContag <- bothInXiProb*log(bothInXiProb)+
                                originOnlyInXi*log(originOnlyInXi)+ 
                                vOnlyInXi*log(vOnlyInXi)+
                                neitherInXi*log(neitherInXi)
  
  return(1/(-log(1/4)) *unnormalisedContag)
}

# my arithmetic tells me it should simplify to
# 1/(-log(1/4)) * (
#   (1+bothInXiProb-2*p)*log(1+bothInXiProb-2*p)
#   + 2*(p-bothInXiProb)*log(p-bothInXiProb)
#   + bothInXiProb*log(bothInXiProb)
# )
#and it matches :)