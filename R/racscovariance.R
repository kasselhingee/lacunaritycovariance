#' @title Covariance Estimation
#' @export racscovariance racscovariance.cvchat
#' @description 
#' Estimates the covariance of a stationary RACS. 
#' The traditional covariance estimator and
#' new estimators based on [**picka2000va**, picka1997va] are available.
#' @author{Kassel Liam Hingee}



#' @param xi A binary map. Either an \code{im} object or a \code{owin} object.
#'   If an \code{im} object then pixel values of 1 or TRUE represent foreground,
#'   0 or \code{FALSE} values represent background, and \code{NA} values
#'   represent outside the  observation window. If an \code{owin} object then
#'   \code{xi} represents foreground and \code{obswin} is required to specify
#'   the observation window.
#' @param obswin The observation window in \code{owin} format if \code{xi} is also in \code{owin} format. 
#' @param setcov_boundarythresh Any vector \eqn{v} such that set covariance of the observation window
#'  is smaller than this threshold is given a covariance of NA to avoid instabilities caused by dividing by very small areas, 
#' @param phat The classical estimate of coverage probability,
#'  which is the observed area in \code{xi} divided by the total area of the observation window.
#'  See \code{coverageprob} for more information.
#' @param cvchat The traditional estimate of covariance in \code{im} format. 
#' Typically created with \code{\link{tradcovarest}}.
#' @param cpp1 Picka's reduced window estimate of coverage probability in \code{im} format - used in improved (balanced) covariance estimators.
#' Can be generated using \code{\link{cppicka}}.
#' @param modification A string specifying the desired modification. See details.
#' @param modifications A list of strings specifying covariance estimators to use. 
#' See details.
#' \code{modifications = "all"} will select all available estimators.  
#' @param drop If TRUE and one modification selected then the returned value will be a single \code{im} object and not a list of \code{im} object.

#' @return 
#' If \code{drop = TRUE} and only one estimator is requested then 
#' an \code{im} object containing the covariance estimate.
#'  Otherwise a named \code{imlist} of covariance estimates corresponding to each requested estimator.


#' @keywords spatial nonparametric

#' @details 
#' The covariance of a RACS is also known as the two-point coverage probability, and is
#' closely related to the semivariogram. 
#' The covariance of a stationary RACS \eqn{\Xi} given a vector \eqn{v} is
#' the probability that two points separated by a vector \eqn{v} are covered by
#' \eqn{\Xi}.
#' 
#' Traditionally covariance for vector \eqn{v} is estimated from a binary map,
#' \eqn{xi}, using the volume of the set of points, \eqn{x}, such that both
#' \eqn{x} and \eqn{x+v} are observed to be in the foreground of \code{xi}
#' relative to the volume of points, \eqn{x}, for which both \eqn{x} and \eqn{x+v}
#' are in the observation window [1].
#' Picka [picka1997va,picka2000va] suggested a number of improvements to centred
#' covariance estimation (see \code{\link{cencovariance}}) that `balanced' the
#' data used to estimate covariance with the data used to estimate coverage
#' probability. These lead to covariance estimators that give
#' estimates for the covariance of \eqn{Xi} that are a constant offset from
#' covariance estimates for the complement of \eqn{Xi} (note the constant offset
#' depends on the coverage probability), which
#' appears to avoid some suprising behaviour that the traditional estimator
#' suffers [Chapter 4, hingee2019thesis].
#' These estimators are called \code{pickaint} and \code{pickaH} in this package.
#' 
#' Another improved estimator inspired by an `intrinsic modification' briefly mentioned by Picka [picka1997va]
#' for pair-correlation estimators`is also available.
#' We have called this estimator \code{mattfeldt} as a similar isotropic estimator for pair-correlation
#' was studied by Mattfeldt and Stoyan [mattfeldt2000im].
#' 

#' The estimators available are (see [Chapter 4, hingee2019thesis] for information): 
#' \itemize{
#' \item{\code{none}} the traditional covariance estimator
#' \item{\code{mattfeldt}} an estimator inspired by an 
#' `intrinsically' balanced pair-correlation estimator from Picka that was later studied in an
#' isotropic situation by Mattfeldt and Stoyan [mattfeldt2000im]
#' \item{\code{pickaint}} an estimator inspired by an 
#' `intrinsically' balanced centred covariance estimator from Picka [picka2000va].
#' \item{\code{pickaH}} an estimator inspired by the 
#' `additively' balanced centred covariance estimator from Picka [picka2000va].
#' }

#' @examples
#' #direct from a binary map
#' xi <- heather$coarse
#' obswin <- Frame(xi)
#' balancedcvchats <- racscovariance(xi, obswin = Frame(xi), modifications = "all")
#' 
#' # from coverage probability estimates and traditional covariance estimate.
#' phat <- coverageprob(xi, obswin = Frame(xi))
#' cvchat <- tradcovarest(xi, inclraw = FALSE)
#' cpp1 <- cppicka(xi, obswin = Frame(heather$coarse))
#' harmonised <- harmonise.im(cvchat = cvchat, cpp1 = cpp1)
#' cvchat <- harmonised$cvchat
#' cpp1 <- harmonised$cpp1
#' 
#' balancedcvchats <- racscovariance.cvchat(cvchat, cpp1, phat, modifications = "pickaH", drop = TRUE)
#' 
#' @describeIn racscovariance Estimates covariance from a binary map.
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

#' @describeIn racscovariance Generates covariances estimates from
#'   a traditional estimate of covariance, Picka's reduced window estimate of coverage probability,
#'   and the traditional estimate of coverage probability.
#'   If these estimates already exist then \code{racscovariance.cvchat} can save significant computation time.
racscovariance.cvchat <- function(cvchat, cpp1 = NULL, phat = NULL, modifications = "all", drop = FALSE){
  harmonised <- harmonise.im(cvchat = cvchat, cpp1 = cpp1)
  cvchat <- harmonised$cvchat
  cpp1 <- harmonised$cpp1
  fcns <- list(
         none = function(cvchat, cpp1 = NULL, phat = NULL) cvchat,
         mattfeldt = balancedracscovariance_mattfeldt_add,
         pickaint = balancedracscovariance_picka_int,
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


balancedracscovariance_mattfeldt_add <- function(cvchat, cpp1, phat){
  return(cvchat - ( (cpp1 + reflect.im(cpp1))/2 )^2 + phat^2) 
}

balancedracscovariance_picka_int <- function(cvchat, cpp1, phat){
  return(cvchat - cpp1*reflect.im(cpp1) + phat^2) 
}

balancedracscovariance_picka_H <- function(cvchat, cpp1, phat){
  return(cvchat - phat*(cpp1 + reflect.im(cpp1) - 2*phat)) 
}

