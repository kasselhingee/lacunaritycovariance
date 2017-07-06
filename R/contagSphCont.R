#' @title Disc State Contagion
#' @export contagSphCont
#' 
#' @description Calculates the disc-state contagion LPI as described in Hingee 2016.
#' The disc-state contagion LPI describes the entropy (mixing) between four possible states of a disc:
#' \enumerate{
#'   \item the disc is completely contained in \eqn{\Xi}
#'   \item the disc does not intersect \eqn{\Xi}
#'   \item the centre of the disc is in \eqn{\Xi} but the disc is not contained in \eqn{Xi}
#'   \item the disc intersects \eqn{Xi} but the centre is outside \eqn{Xi}
#' }
#' The user is required to specify the radius of the disc which is essentially a distance for which the user is interested in quantifying mixing. 
#' 
#' The difference to classical contagion is that disc-state contagion is based on the spherical contact distribution instead of pixel neighbours.
#' One impact of this design is that a distance to quantify the mixing between \eqn{Xi} and the background is chosen by the user (for classical contagion this distance is fixed by the image resolution).
#' 
#' @param XiH Conditional spherical contact distribution function for \eqn{\Xi}. 
#' Typically this is an \code{fv} object but could also be a vector of values.
#' In applications \code{XiH} would likely be estimated from an image using \code{\link{Hest}} in \pkg{spatstat}.
#' @param XiHc Conditional spherical contact distribution for the complement of \eqn{\Xi}. 
#' This is called the Conditional Core Probability in Hingee 2016.
#' Typically this is an \code{fv} object but could also be a vector of values.
#' In applications \code{XiH} would likely be estimated from an image using \code{\link{Hest}} in \pkg{spatstat}.
#' @param p  The coverage fraction of \eqn{\Xi}.
#' In applications to images an estimate of the coverage fraction can be obtained using \code{\link{coveragefrac}}.
#' @param normalise Optional. If TRUE \code{contagSphCont} normalises the results so that all RACS return a value between 0 and 1. Default is FALSE. 
#' @details XiH should be a function of radius that gives (or estimates) the probability of a disc of radius \eqn{r} not intersecting \eqn{\Xi} if the disc's centre is not in \eqn{\Xi} 
#' \deqn{XiH(r) = P(B_r(x) \subseteq \Xi^c | x \in \Xi^c).}
#' Similarly \code{XiHc} should be an estimate of the probability of a disc being fully contained in \eqn{\Xi} given its centre is in \eqn{\Xi}
#' \deqn{XiHc(r)\approx P(B_r(x) \subseteq \Xi | x \in \Xi).}
#' These can both be obtained using \code{\link{Hest}} in \pkg{spatstat}.
#' For \code{XiHc} take care to apply Hest to the complement of \eqn{\Xi} with the observation window \eqn{W}.
#' 
#' If \code{XiH} and \code{XiHc} are both fv objects then they must be generated using Hest because the function automatically uses the reduce-sample border correction estimates.
#' In this case the return value is an fv object.
#'
#' If \code{normalise} is \code{TRUE} then the result is divided by 
#' \eqn{\frac{-4}{e}ln(\frac{1}{e})} and added to 1 so that the normalised disc state contagion is between 0 and 1.
#'
#' @return An \code{fv} object or a vector the same length as \code{XiH} corresponding to the contagion at each r value of \code{XiH}

#' @references 
#' Hingee, K.L. (2016) Statistics for Patch Observations. ISPRS Congress Proceedings p. IPSRS.


#' @examples
#' xi <- heather$coarse
#' obswindow <- Frame(heather$coarse)
#' p <- coveragefrac(xi,Frame(xi))
#' XiH <- Hest(xi,W=obswindow) #Sph. Contact Distrution Estimate
#' XiHc <- Hest(complement.owin(xi),W=obswindow) #Conditional Core Prob. Estimate
#' plot(XiH,type="l",col="red") 
#' lines(XiHc,type="l",col="black") 
#' 
#' contagion <- contagSphCont(XiH,XiHc,p)
#' plot(contagion)
#' 

contagSphCont <- function(XiH, XiHc, p, normalise=FALSE){
  returnfv <- FALSE
  unitnames <- NULL
  if (is.fv(XiH) && is.fv(XiHc)) {#then new version of contagion
    returnfv <- TRUE
    unitnames <- unitname(XiH)
    fvin <- list(XiH=XiH,XiHc=XiHc)
    XiHf <- as.function.fv(XiH,value=".y",extrapolate=TRUE)
    XiHcf <- as.function.fv(XiHc,value=".y",extrapolate=TRUE)
    argranges <- lapply(list(XiH=XiH,XiHc=XiHc),argumentrange)
    ## determine finest resolution (from spatstat)
    xsteps <- sapply(fvin, argumentstep)
    finest <- which.min(xsteps)
    ## extract argument values (from spatstat)
    xvals <- with(fvin[[finest]], .x)
    fvwlargestarg <- which.max(list(XiH=argranges$XiH[[2]],XiHc=argranges$XiHc[[2]]))
    xvals <- c(xvals,with(fvin[[fvwlargestarg]], .x)[with(fvin[[fvwlargestarg]], .x)>max(xvals)])
    ##

    XiH <- XiHf(xvals)
    XiHc <- XiHcf(xvals)
  }
  Pstates <- matrix(NA,nrow=4,ncol=length(XiH))
  rownames(Pstates)=c("P11","P10","P01","P00")
  Pstates["P11",] <- p* (1-XiHc)
  Pstates["P10",] <- p*XiHc
  Pstates["P00",] <- (1-p)*(1-XiH)
  Pstates["P01",] <- (1-p)*(XiH)
  
  tempstates <- Pstates
  tempstates[Pstates<1E-8] <- 1
  contag <- colSums(Pstates*log(tempstates))
  if (normalise) {contag <- 1+ contag/(-4/exp(1)*log(1/exp(1)))}
  if (returnfv){
    if (normalise) {
      return(fv(data.frame(r= xvals,
                           contag = contag),
                valu = "contag",
                desc=c("radius",
                       "normalised SCD contagion estimate"),
                unitname=unitnames
      ))
    }
    else {
      return(fv(data.frame(r= xvals,
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


argumentrange <- function(f) { range(with(f, .x)) }  #copied from spatstat code
argumentstep <- function(f) { mean(diff(with(f, .x))) }  #copied from spatstat code
