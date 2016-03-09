#' @title Spherical Contact version of Contagion
#' @export contagSphCont
#' 
#' @description Calculates a spherical contact version of contagion. 
#' This version of contagion gauges closeness/mixing between \eqn{\Xi} and \eqn{\Xi^c}
#'  using the spherical contact distribution. 
#'  Given a distance \eqn{r} and an arbitrary point \eqn{x}
#'  in \eqn{\Xi} then \eqn{\Xi} close means everything with \eqn{r} of \eqn{x} is in \eqn{\Xi}.
#'  \eqn{\Xi^c} close to \eqn{x} means there is at least some part of \eqn{\Xi^c} within \eqn{r} of \eqn{x}.
#' 
#' @param xiH Spherical contact function for xi. 
#' @param xiHc Spherical contact distribution for the complement of xi.
#' If xiH and xiHc are both fv objects then harmonises them and function returns an fv object
#' @param p  An estimate of the coverage fraction of a RACS \eqn{\Xi}.
#'  The sample points should be same as xiH for now.
#' @param normalise Optional. If true normalises the results so that all RACS return a value between 
#' @details xiH should be a function of radius that estimates the probability 
#' \deqn{xiH(r)\approx P(B_r(x) \subseteq \Xi^c)}
#'  of a disc around an arbitrary point \eqn{x} is contained in \eqn{\Xi^c.}
#' Similary xiHc should be the probability of disc being fully contained in \eqn{\Xi}
#' \deqn{xiHc(r)\sim P(B_r(o) \subseteq \Xi).}
#' 
#' If \code{normalise} is \code{TRUE} then divides by 
#' \eqn{\frac{-4}{e}ln(\frac{1}{e})} and adds 1 so normalised spherical contact contagion is
#' \deqn{
#' 1+(\frac{-4}{e}ln(\frac{1}{e}))^{-1} \mbox{unnormalised contagion}
#' }
#' This makes contagion vary between 0 and 1 for all 2 phase processes.
#' @return a vector the same length as xiH corresponding to the contagion at each r value of xiH


#' @examples
#' xi <- heather$coarse
#' p <- coveragefrac(xi,Frame(xi))
#' xiH <- Hest(xi)
#' xiHc <- Hest(complement.owin(xi))
#' #it typically not advisable to choose set r ourselves, **is it a good idea here? Interpolation later might be better?
#' plot(xiH,type="l",col="red") 
#' lines(xiHc,type="l",col="black") 
#' 
#' harmonised <- harmonise(xiH,xiHc)
#' xiH <- harmonised[[1]]
#' xiHc <- harmonised[[2]]
#' 
#' contagion <- contagSphCont(xiH$km,xiHc$km,p)
#' plot(xiH$r,contagion,type="l")
#' 
#' contagion <- contagSphCont(xiH,xiHc,p)
#' plot(contagion)
#' 
#' @seealso \code{\link{contagTwoPtProb}} 

contagSphCont <- function(xiH, xiHc, p, normalise=FALSE){
  returnfv <- FALSE
  if (is.fv(xiH) && is.fv(xiHc)) {#then new version of contagion
    returnfv <- TRUE
    harmonisedSCDs <- harmonise(xiH,xiHc)
    r <- harmonisedSCDs[[1]]$r
    xiH <- harmonisedSCDs[[1]]$rs
    xiHc <- harmonisedSCDs[[2]]$rs
  }
  Pstates <- matrix(NA,nrow=4,ncol=length(xiH))
  rownames(Pstates)=c("P11","P10","P01","P00")
  Pstates["P11",] <- p* (1-xiHc)
  Pstates["P10",] <- p*xiHc
  Pstates["P00",] <- (1-p)*(1-xiH)
  Pstates["P01",] <- (1-p)*(xiH)
  
  tempstates <- Pstates
  tempstates[Pstates<1E-8] <- 1
  contag <- colSums(Pstates*log(tempstates))
  if (normalise) {contag <- 1+ unnormalisedContag/(-4/exp(1)*log(1/exp(1)))}
  if (returnfv){
    if (normalise) {
      return(fv(data.frame(r= r,
                           contag = contag),
                valu = "contag",
                desc=c("radius",
                       "normalised SCD contagion estimate")
      ))
    }
    else {
      return(fv(data.frame(r= r,
                           contag = contag),
                valu = "contag",
                desc=c("radius",
                       "unnormalised SCD contagion estimate")
      ))
    }
  }
  return(contag)
}

