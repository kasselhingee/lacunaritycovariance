#' @title Spherical Contact version of Contagion
#' @export contagSphCont
#' 
#' @param xiH Spherical contact function for xi. 
#' @param xiHc Spherical contact distribution for the complement of xi.
#' @param p  An estimate of the coverage fraction of a RACS $\Xi$.
#'  The sample points should be same as xiH for now.
#' @param normalise Optional. If true normalises the results so that all RACS return a value between 
#' @details xiH should be a function of radius that estimates the probability 
#' \deqn{xiH(r)\approx P(B_r(x) \subseteq \Xi^c)}
#'  of a disc around an arbitrary point $x$ is contained in \eqn{\Xi^c.}
#' Similary xiHc should be the probability of disc being fully contained in $\Xi$
#' \deqn{xiHc(r)\approx P(B_r(o) \subseteq \Xi).}
#' 
contagSphCont <- function(xiH,xiHc,p=NULL,normalise=FALSE){
  Pstates <- matrix(NA,nrow=4,ncol=length(xiH))
  rownames(Pstates)=c("P11","P10","P01","P00")
  Pstates["P11",] <- p* (1-xiHc)
  Pstates["P10",] <- p*xiHc
  Pstates["P00",] <- (1-p)*(1-xiH)
  Pstates["P01",] <- (1-p)*(xiH)
  
  tempstates <- Pstates
  tempstates[Pstates<1E-8] <- 1
  unnormalisedContag <- colSums(Pstates*log(tempstates))
  if (normalise) {return(1+ (4/exp(1))*unnormalisedContag)}
  else {return(unnormalisedContag)}
}

#' @examples 
#' xi <- heather$coarse
#' p <- covpest(xi,Frame(xi))
#' xiH <- Hest(xi)
#' r <- xiH$r
#' xiH <- xiH$km
#' #it typically not advisable to choose set r ourselves, **is it a good idea here? Interpolation later might be better?
#' xiHc <- Hest(complement.owin(xi,frame=Frame(xi)),r=r)$km
#' plot(r,xiH,type="l",col="red") 
#' lines(r,xiHc,type="l",col="black") 
#' 
#' contagion <- contagSphCont(xiH,xiHc,p)
#' plot(r,contagion,type="l")
#' 