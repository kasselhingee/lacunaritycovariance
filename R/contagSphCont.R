#' @title Disc State Contagion
#' @export scdcontagion
#' 
#' @description Calculates the disc-state contagion LPI as described in Hingee 2016.
#' The disc-state contagion LPI describes the entropy (mixing) between four possible states of a disc:
#' \enumerate{
#'   \item the disc is completely contained in \eqn{\Xi}
#'   \item the disc does not intersect \eqn{\Xi}
#'   \item the centre of the disc is in \eqn{\Xi} but the disc is not contained in \eqn{\Xi}
#'   \item the disc intersects \eqn{\Xi} but the centre is outside \eqn{\Xi}
#' }
#' 
#' Dics-state contagion is a function of the disc radius.
#' 
#' The main difference to classical contagion is that disc-state contagion is based on the spherical contact distribution instead of pixel neighbours.
#' One impact of this design is that the distance with which to quantify the mixing between \eqn{\Xi} and the background may be chosen by the user by choosing the disc radius (for classical contagion this distance is fixed by the image resolution).
#' 
#' @param Xi.H Conditional spherical contact distribution function for \eqn{\Xi}. 
#' Typically this is an \code{fv} object but could also be a vector of values.
#' In applications \code{Xi.H} would likely be estimated from an image using \code{\link{Hest}} in \pkg{spatstat}.
#' @param Xic.H Conditional spherical contact distribution for the complement of \eqn{\Xi}. 
#' This is called the Conditional Core Probability in Hingee 2016.
#' Typically this is an \code{fv} object but could also be a vector of values.
#' In applications \code{Xi.H} would likely be estimated from an image using \code{\link{Hest}} in \pkg{spatstat}.
#' @param p  The coverage probability of \eqn{\Xi}.
#' In applications to images an estimate of the coverage probability can be obtained using \code{\link{coverageprob}}.
#' @param normalise Optional. If TRUE \code{scdcontagion} normalises the results so that all RACS return a value between 0 and 1. Default is FALSE. 
#' @details Xi.H should be a function of radius that gives (or estimates) the probability of a disc of radius \eqn{r} not intersecting \eqn{\Xi} if the disc's centre is not in \eqn{\Xi} 
#' \deqn{\code{Xi.H}(r) = P(B_r(x) \subseteq \Xi^c | x \in \Xi^c).}
#' Similarly \code{Xic.H} should be an estimate of the probability of a disc being fully contained in \eqn{\Xi} given its centre is in \eqn{\Xi}
#' \deqn{\code{Xic.H}(r)\approx P(B_r(x) \subseteq \Xi | x \in \Xi).}
#' These can both be obtained using \code{\link{Hest}} in \pkg{spatstat}.
#' For \code{Xic.H} take care to apply Hest to the complement of \eqn{\Xi} with the observation window \eqn{W}.
#' 
#' If \code{normalise} is \code{TRUE} then the result is divided by 
#' \eqn{\frac{-4}{e}ln(\frac{1}{e})} and added to 1 so that the return value is between 0 and 1 for all possible \eqn{\Xi}.
#'
#' @return An \code{fv} object or a vector the same length as \code{Xi.H} corresponding to the contagion at each r value of \code{Xi.H}

#' @references 
#' Hingee, K.L. (2016) Statistics for Patch Observations. ISPRS - International Archives of the Photogrammetry, Remote Sensing and Spatial Information Sciences pp. 235-242. ISPRS.


#' @examples
#' xi <- heather$coarse
#' obswindow <- Frame(heather$coarse)
#' p <- coverageprob(xi,Frame(xi))
#' xi.H <- Hest(xi,W=obswindow) #Sph. Contact Distrution Estimate
#' xic.H <- Hest(complement.owin(xi),W=obswindow) #Conditional Core Prob. Estimate
#' plot(xi.H,type="l",col="red") 
#' lines(xic.H,type="l",col="black") 
#' 
#' contagion <- scdcontagion(xI.H,xic.H,p)
#' plot(contagion)
#' 
#' @keywords spatial nonparametric 
scdcontagion <- function(Xi.H, Xic.H, p, normalise=FALSE){
  returnfv <- FALSE
  unitnames <- NULL
  if (is.fv(Xi.H) && is.fv(Xic.H)) {
    returnfv <- TRUE
    unitnames <- unitname(Xi.H)
    fvin <- list(Xi.H = Xi.H, Xic.H = Xic.H)
    Xi.Hf <- as.function.fv(Xi.H, value = ".y", extrapolate = TRUE)
    Xic.Hf <- as.function.fv(Xic.H, value = ".y", extrapolate = TRUE)
    argranges <- lapply(list(Xi.H = Xi.H, Xic.H = Xic.H), argumentrange)
    ## determine finest resolution (from spatstat)
    xsteps <- sapply(fvin, argumentstep)
    finest <- which.min(xsteps)
    ## extract argument values (from spatstat)
    xvals <- with(fvin[[finest]], .x)
    fvwlargestarg <- which.max(list(Xi.H = argranges$Xi.H[[2]],
                                   Xic.H = argranges$Xic.H[[2]]))
    xvals <- c(xvals,
      with(fvin[[fvwlargestarg]], .x)[
        with(fvin[[fvwlargestarg]], .x) > max(xvals)
      ])
    ##

    Xi.H <- Xi.Hf(xvals)
    Xic.H <- Xic.Hf(xvals)
  }
  probofstate <- matrix(NA, nrow = 4, ncol = length(Xi.H))
  rownames(probofstate) <- c("P11", "P10", "P01", "P00")
  probofstate["P11", ] <- p * (1 - Xic.H)
  probofstate["P10", ] <- p * Xic.H
  probofstate["P00", ] <- (1 - p) * (1 - Xi.H)
  probofstate["P01", ] <- (1 - p) * (Xi.H)

  tempstates <- probofstate
  tempstates[probofstate < 1E-8] <- 1
  contag <- colSums(probofstate * log(tempstates))
  if (normalise) {
    contag <- 1 + contag / (-4 / exp(1) * log(1 / exp(1)))
  }
  if (returnfv) {
    if (normalise) {
      return(fv(data.frame(r = xvals,
                           contag = contag),
                valu = "contag",
                desc = c("radius",
                       "normalised SCD contagion estimate"),
                unitname = unitnames
      ))
    }
    else {
      return(fv(data.frame(r = xvals,
                           contag = contag),
                valu = "contag",
                desc = c("radius",
                       "unnormalised SCD contagion estimate"),
                unitname = unitnames
      ))
    }
  }
  return(contag)
}

#copied from spatstat code
argumentrange <- function(f) {
  range(with(f, .x))
}
#copied from spatstat code
argumentstep <- function(f) {
  mean(diff(with(f, .x)))
}
