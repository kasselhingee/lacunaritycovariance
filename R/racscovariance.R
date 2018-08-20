#' @title Covariance Estimation
#' @export racscovariance balancedracscovariance.cvchat  racscovariance.cvchat
#' @description 
#' Estimates the covariance of a stationary RACS from a binary map. 
#' The traditional covariance estimator and
#' new estimators based on [**picka2000va**, picka1997va] are available.
#' @author{Kassel Liam Hingee}



#' @param xi A binary map. Either an \code{im} object with pixel values of 0, 1, TRUE, FALSE or NA, 
#' or an \code{owin} object representing the foreground (in which case \code{obswin} is required).
#' @param obswin The observation window in \code{owin} format if \code{xi} is also in \code{owin} format. 
#' @param setcov_boundarythresh Any vector \eqn{v} such that set covariance of the observation window
#'  is smaller than this threshold is given a covariance of NA to avoid instabilities caused by dividing by very small areas, 
#' @param phat The classical estimate of coverage probability,
#'  which is the observed area in \code{xi} divided by the total area of the observation window.
#'  See \code{coverageprob} for more information.
#' @param cvchat The traditional estimate of covariance in \code{im} format. 
#' Typically created with \code{\link{tradcovarest}}.
#' @param cpp1 Picka's reduced estimate of coverage probability in \code{im} format - used in improved (balanced) covariance estimators.
#' Can be generated using \code{\link{cppicka}}.
#' @param modification A string specifying the desired modification. See details.
#' @param modifications A list of strings specifying covariance estimators to use. 
#' See details.
#' \code{modifications = "all"} will select all available estimators.  
#' @param drop If TRUE and one modification selected then the returned value will be a single \code{im} object and not a list of \code{im} object.

#' @return 
#' If \code{drop = TRUE} and only one estimator is requested then 
#' an \code{im} object containing the covariance estimate.
#'  Otherwise a named \code{imlist} of such covariance estimates corresponding to each requested estimator.


#' @keywords spatial nonparametric

#' @details 
#' The covariance of a RACS is also known as the two-point coverage probability, and is
#' closely related to the semivariogram. 
#' The covariance of a stationary RACS \eqn{\Xi} given a vector \eqn{v} is
#' the probability that two points separated by a vector \eqn{v} are covered by
#' \eqn{\Xi}.
#' 
#' Traditionally covariance for vector \eqn{v} is estimated from a binary map
#' \eqn{xi} using the volume of the set of points, \eqn{x}, such that both
#' \eqn{x} and \eqn{x+v} are observed to be in the foreground of \code{xi}
#' relative to the volume of points, \eqn{x}, for which both \eqn{x} and \eqn{x+v}
#' are in the observation window [1].
#' 
#'   ,
#' \deqn{C_T (v) = \gamma_{\cap (W, X)}(v)/ \le \gamma_W(v)}
#'  where
#' \eqn{X} is a realisation of \eqn{\Xi} observed in window \eqn{W},
#' \eqn{\gamma_{W}(v)} is the set covariance of the observation window \eqn{|W
#' \cap (W\oplus v)|} and \eqn{\gamma_{W\cap X}(v)} is the set covariance of an
#' observation, \eqn{X}, of the RACS \eqn{\Xi}, \eqn{\gamma_{W\cap X}(v) = |W
#' \cap X \cap ((W\cap X ) \oplus v)|}.
#' 
#' Picka [picka1997va,picka2000va] suggested a number of improvements to
#'  centred covariance estimation
#'   (see \code{\link{cencovariance}}) that
#'    `balanced' the use of data in the observation window.
#'  These lead to estimators \eqn{\hat{C}(v)} of covariance that satisify
#'  \deqn{\hat{C}(v, X) = \hat{C}(v, X^c) - 1 + 2\hat{p}},
#'  which is the same as the relationship between 
#'  covariance of \eqn{\Xi}, \eqn{\Xi^c} and the coverage probability of \eqn{\Xi}.
#'  It appears that these new estimators of covariance avoid
#'   some suprising behaviour that the traditional estimator has [Chapter 4, hingee2019thesis].
#' 
#' If \code{xi} is an \code{im} object the pixel values of 1 or TRUE are interpreted as foreground
#'  and the pixel values of 0 or FALSE are interpreted as background.
#'  Pixel values of NA are considered outside the observation window.
#' 
#' Modifies the classical covariance estimator to based on balancing ideas for pair-correlation and centred covariance estimation.
#' Many of the modifications use Picka's coverage probability estimators in some way.
#' Modifications available are: 
#' \itemize{
#' \item \code{none} Returns cvchat
#' \item symm Returns the estimated average of cvchat at -v and +v
#' \item adrian  Similar to Picka additive but uses \code{cpp1} twice and not the refelction of \code{cpp1}
#' \item mattfeldadd
#' \item mattfeldmult
#' \item pickaint The only modification supplied that provides estimates that satisfy C(v) = 2phat - 1 + Ccomplment(v)
#' \item pickaintmult
#' \item pickaH
#' }

#' @examples
#' xi <- heather$coarse
#' obswin <- Frame(xi)
#' balancedcvchats <- racscovariance(xi, obswin = Frame(xi), modifications = "all")
#' balancedcvchats2 <- byconv_cvchats(xi, obswin = Frame(xi), modifications = "all")
#' 
#' # plot.solist(c(balancedcvchats, balancedcvchats2), ncols = 8)
#' diff <- mapply(function(x, y) x - y, balancedcvchats, balancedcvchats2, SIMPLIFY = FALSE)
#' # plot.solist(diff)
#' 

#' phat <- coverageprob(xi, obswin = Frame(xi))
#' cvchat <- tradcovarest(xi, inclraw = FALSE)
#' cpp1 <- cppicka(xi, obswin = Frame(heather$coarse))
#' harmonised <- harmonise.im(cvchat = cvchat, cpp1 = cpp1)
#' cvchat <- harmonised$cvchat
#' cpp1 <- harmonised$cpp1
#' 
#' balancedcvchat <- balancedracscovariance.cvchat(cvchat, cpp1, phat, modification = "pickaint")
#' balancedcvchats <- racscovariance.cvchat(cvchat, cpp1, phat, modifications = "all")
#' modifications <- c("none",
#'  "symm",
#'  "adrian",
#'  "mattfeldt",
#'  "mattfeldtmult",
#'  "pickaint",
#'  "pickaintmult",
#'  "pickaH", 
#'  function(cvchat, cpp1, phat) cvchat)
#' balancedcvchats <- racscovariance.cvchat(cvchat, cpp1, phat, modifications = modifications)
#' # plot(as.solist(balancedcvchats), equal.ribbon = TRUE)
#' 

#' 
racscovariance <- function(xi, obswin = NULL,
        setcov_boundarythresh = NULL,
        modifications = "all",
        drop = FALSE){
  cvchat <- tradcovarest(xi, obswin, setcov_boundarythresh = setcov_boundarythresh)
  cpp1 <- cppicka(xi, obswin, setcov_boundarythresh = setcov_boundarythresh)
  phat <- coverageprob(xi, obswin)
  
  cvchats <- racscovariance.cvchat(cvchat, cpp1, phat, modifications = modifications, drop = drop) 
  return(cvchats)
}


#' @describeIn racscovariance Applies covariance balancing modification to precomputed cvchat, cpp1 and phat
balancedracscovariance.cvchat <- function(cvchat, cpp1 = NULL, phat = NULL, modification = "pickaH"){
  harmonised <- harmonise.im(cvchat = cvchat, cpp1 = cpp1)
  cvchat <- harmonised$cvchat
  cpp1 <- harmonised$cpp1
  balancedcvchat <- switch(modification,
         none = cvchat,
         symm = balancedracscovariance_symm(cvchat),
         adrian = balancedracscovariance_adrian(cvchat, cpp1, phat),
         mattfeldt = balancedracscovariance_mattfeldt_add(cvchat, cpp1, phat),
         mattfeldtmult = balancedracscovariance_mattfeldt_mult(cvchat, cpp1, phat),
         pickaint = balancedracscovariance_picka_int(cvchat, cpp1, phat),
         pickaintmult = balancedracscovariance_picka_intmult(cvchat, cpp1, phat),
         pickaH = balancedracscovariance_picka_H(cvchat, cpp1, phat),
         stop(paste("Modification", modification, "not found.")) 
         )
  return(balancedcvchat)
}

#' @describeIn racscovariance Applies multiple modifications simultaneously from a precomputed cvchat, cpp1 and phat
racscovariance.cvchat <- function(cvchat, cpp1 = NULL, phat = NULL, modifications = "all", drop = FALSE){
  harmonised <- harmonise.im(cvchat = cvchat, cpp1 = cpp1)
  cvchat <- harmonised$cvchat
  cpp1 <- harmonised$cpp1
  fcns <- list(
         none = function(cvchat, cpp1 = NULL, phat = NULL) cvchat,
         symm = balancedracscovariance_symm,
         adrian = balancedracscovariance_adrian,
         mattfeldt = balancedracscovariance_mattfeldt_add,
         mattfeldtmult = balancedracscovariance_mattfeldt_mult,
         pickaint = balancedracscovariance_picka_int,
         pickaintmult = balancedracscovariance_picka_intmult,
         pickaH = balancedracscovariance_picka_H
  )
  if ((modifications == "all")[[1]]) {modifications <- names(fcns)}
  fcnstouse <- fcns[names(fcns) %in% modifications]
  isfunction <- unlist(lapply(modifications, function(x) "function" %in% class(x)))
  modificationsnotused <- modifications[!( (modifications %in% names(fcns)) | isfunction)]
  
  fcnstouse <- c(fcnstouse, modifications[isfunction]) #add user specified modification
  
  if(length(modificationsnotused) > 0){stop(
    paste("The following modifications are not recognised as existing function names or as a function:", modificationsnotused))}
  balancedcvchats <- lapply(fcnstouse, function(x) do.call(x, args = list(cvchat = cvchat, cpp1 = cpp1, phat = phat)))

  if (drop && (length(balancedcvchats) == 1)){ return(balancedcvchats[[1]])
  } else {return(as.imlist(balancedcvchats)) }
}


balancedracscovariance_symm <- function(cvchat, cpp1 = NULL, phat = NULL){
  return((cvchat + reflect.im(cvchat))/2) 
}

balancedracscovariance_adrian <- function(cvchat, cpp1, phat){
  return(cvchat - cpp1*cpp1 + phat^2) 
}

balancedracscovariance_mattfeldt_add <- function(cvchat, cpp1, phat){
  return(cvchat - ( (cpp1 + reflect.im(cpp1))/2 )^2 + phat^2) 
}

balancedracscovariance_mattfeldt_mult <- function(cvchat, cpp1, phat){
  return(cvchat * phat^2/ ( ( (cpp1 + reflect.im(cpp1))/2 )^2) ) 
}

balancedracscovariance_picka_int <- function(cvchat, cpp1, phat){
  return(cvchat - cpp1*reflect.im(cpp1) + phat^2) 
}

balancedracscovariance_picka_intmult <- function(cvchat, cpp1, phat){
  return(cvchat * phat^2 / (cpp1*reflect.im(cpp1))) 
}

balancedracscovariance_picka_H <- function(cvchat, cpp1, phat){
  return(cvchat - phat*(cpp1 + reflect.im(cpp1) - 2*phat)) 
}

