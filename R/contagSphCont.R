#' @title Disc State Contagion
#' @export contagSphCont
#' 
#' @description Calculates the disc-state contagion as described in Hingee 2016. It is like the contagion LPI but is based on the spherical contact version of contagion. 
#' It describes the entropy (mixing) between four possible states of disc
#' (see Hingee 2016 for more details).
#' It requires a mixing distance of interest to be chosen by the user
#' (compared to classical contagion for which this distance is set by the image resolution).
#' 
#' @param xiH Estimated conditional spherical contact distribution function for \eqn{\Xi}. 
#' Typically as a \code{fv} object but could also be a vector of values.
#' @param xiHc Estimate conditional spherical contact distribution for the complement of \eqn{\Xi}. 
#' This is called the Conditional Core Probability in Hingee 2016.
#' Typically is a \code{fv} object.
#' @param p  An estimate of the coverage fraction of a RACS \eqn{\Xi}.
#' Typically obtained using \code{coveragefrac}.
#' @param normalise Optional. If TRUE \code{contagSphCont} normalises the results so that all RACS return a value between 0 and 1. Default is FALSE. 
#' @details xiH should be a function of radius that estimates the probability 
#' \deqn{xiH(r)\approx P(B_r(x) \subseteq \Xi^c | x \in \Xi^c)}
#'  of a disc around an arbitrary point \eqn{x} is contained in \eqn{\Xi^c.}
#' Similary xiHc should be an estimate of the probability of a disc being fully contained in \eqn{\Xi}
#' \deqn{xiHc(r)\approx P(B_r(x) \subseteq \Xi | x \in \Xi).}
#' These can both be obtained using \code{Hest} in \code{spatstat}.
#' 
#' If xiH and xiHc are both fv objects then they must be generated using Hest because the function automatically uses the reduce-sample border correction estimates.
#' In this case the return value is an fv object.
#'
#' If \code{normalise} is \code{TRUE} then the result is divided by 
#' \eqn{\frac{-4}{e}ln(\frac{1}{e})} and added to 1 so that the normalised disc state contagion is between 0 and 1.
#'
#' @return An \code{fv} object or a vector the same length as xiH corresponding to the contagion at each r value of xiH

#' @references 
#' Hingee, K.L. (2016) Statistics for Patch Observations. ISPRS Congress Proceedings p. IPSRS.


#' @examples
#' xi <- heather$coarse
#' obswindow <- Frame(heather$coarse)
#' p <- coveragefrac(xi,Frame(xi))
#' xiH <- Hest(xi,W=obswindow) #Sph. Contact Distrution Estimate
#' xiHc <- Hest(complement.owin(xi),W=obswindow) #Conditional Core Prob. Estimate
#' plot(xiH,type="l",col="red") 
#' lines(xiHc,type="l",col="black") 
#' 
#' contagion <- contagSphCont(xiH,xiHc,p)
#' plot(contagion)
#' 

contagSphCont <- function(xiH, xiHc, p, normalise=FALSE){
  returnfv <- FALSE
  unitnames <- NULL
  if (is.fv(xiH) && is.fv(xiHc)) {#then new version of contagion
    returnfv <- TRUE
    unitnames <- unitname(xiH)
    harmonisedSCDs <- harmonise(xiH,xiHc)
    r <- harmonisedSCDs[[1]]$r
    rharmleng <- length(r)
    #the following if/else statements extends the harmonised values when its known that xiH==1 or xiHc==1
    if (max(xiH$r)>r[rharmleng] & (harmonisedSCDs[[2]]$rs[rharmleng]>0.99)){
      r <- c(r,xiH$r[xiH$r>r[rharmleng]])
      xiH <- c(harmonisedSCDs[[1]]$rs,xiH$rs[xiH$r>r[rharmleng]])
      xiHc <- c(harmonisedSCDs[[2]]$rs,rep(1,length(r)-rharmleng))
    }
    else if (max(xiHc$r)>r[rharmleng] & (harmonisedSCDs[[1]]$rs[rharmleng]>0.99)){
      r <- c(r,xiHc$r[xiHc$r>r[rharmleng]])
      xiHc <- c(harmonisedSCDs[[1]]$rs,xiHc$rs[xiHc$r>r[rharmleng]])
      xiH <- c(harmonisedSCDs[[2]]$rs,rep(1,length(r)-rharmleng))
    }
    else {
      xiH <- harmonisedSCDs[[1]]$rs
      xiHc <- harmonisedSCDs[[2]]$rs
    }
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
  if (normalise) {contag <- 1+ contag/(-4/exp(1)*log(1/exp(1)))}
  if (returnfv){
    if (normalise) {
      return(fv(data.frame(r= r,
                           contag = contag),
                valu = "contag",
                desc=c("radius",
                       "normalised SCD contagion estimate"),
                unitname=unitnames
      ))
    }
    else {
      return(fv(data.frame(r= r,
                           contag = contag),
                valu = "contag",
                desc=c("radius",
                       "unnormalised SCD contagion estimate"),
                unitname=unitnames
      ))
    }
  }
  return(contag)
}

