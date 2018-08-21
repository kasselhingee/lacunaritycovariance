#' @title Balanced estimation of pair-correlation.
#' @export paircorr  paircorr.cvchat
#' @description 
#' Estimates the pair-correlation function of a stationary RACS. 
#' The traditional pair-correlation estimator and three 'balanced' estimators suggested by Picka
#' are available.
#' @author{Kassel Liam Hingee}

#' @param xi A binary map of an observation of a RACS of interest. See
#'   \code{\link{stationaryracsinference-package}} for details.
#' @param obswin If \code{xi} is an \code{owin} object then \code{obswin} is an
#'   \code{owin} object that specifies the observation window.
#' @param setcov_boundarythresh Any vector \eqn{v} such that set covariance of the observation window
#'  is smaller than this threshold is given a covariance of NA to avoid instabilities caused by dividing by very small areas, 
#' @param phat The traditional estimate of coverage probability,
#'  which is the observed foreground area in \code{xi} divided by the total area of the observation window.
#'  See \code{\link{coverageprob}} for more information.
#' @param cvchat The traditional estimate of covariance in \code{im} format. 
#' Typically created with \code{\link{tradcovarest}}.
#' @param cpp1 Picka's reduced window estimate of coverage probability in \code{im} format - used in improved (balanced) covariance estimators.
#' Can be generated using \code{\link{cppicka}}.
#' @param modifications A list of strings specifying estimators to use. 
#' See details.
#' \code{modifications = "all"} will select all available estimators.  
#' @param drop If TRUE and one modification selected then the returned value will be a single \code{im} object and not a list of \code{im} object.
#'  modifications = "all" will select all inbuilt modifications. See details. 

#' @return If \code{drop = TRUE} and a single estimator requested then a
#'   \code{im} object containing the pair-correlaion estimate. Otherwise a
#'   named \code{imlist} of \code{im} objects containing the pair-correlation
#'   estimates for each requested estimator.
#'
#'   


#' @keywords spatial nonparametric
#' @details The pair-correlation of a stationary RACS is 
#' \deqn{g(v) = C(v) / p^2.}
#'
#'   The estimators available are (see [Chapter 4, hingee2019thesis] for
#'   more information): 
#'   \itemize{ 
#'   \item{\code{none}} the traditional pair-correlation estimator which is \eqn{Chat(v) / (phat^2)}, where \eqn{Chat} and \eqn{phat} are 
#' the traditional estimates of coverage probability and covariance respectively. 
#'   \item{\code{mattfeldt}} an `intrinsically' balanced pair-correlation estimator suggested by Picka.
#'   A similar isotropic pair-correlation estimator was later studied by Mattfeldt and Stoyan [mattfeldt2000im].
#'   \item{\code{pickaint}} Picka's 'intrinsically' balanced pair-correlation estimator [picka2000va]. 
#'   \item{\code{pickaH}} Picka's 'additively' balanced pair-correlation estimator [picka2000va].
#'   }
#'
#' @examples
#' xi <- heather$coarse
#' #estimate directly from a binary map
#' pclns_direst <- paircorr(as.im(xi, na.replace = 0), modifications = "all")
#' 
#' #estimate using traditional covariance estimates, traditional coverage
#' probability estimate and Picka's reduced window coverage probability estimates.
#' obswin <- Frame(xi)
#' phat <- coverageprob(xi, obswin = Frame(xi))
#' cvchat <- tradcovarest(xi)
#' cpp1 <- cppicka(xi, obswin = Frame(heather$coarse))
#' pclns_frcvc <- paircorr.cvchat(cvchat, cpp1, phat, modifications = "all")

#' @describeIn paircorr Estimates pair-correlation from a binary map.
paircorr <- function(xi, obswin = NULL,
                  setcov_boundarythresh = NULL,
                  modifications = "all",
                  drop = FALSE){
  cvchat <- tradcovarest(xi, obswin, setcov_boundarythresh = setcov_boundarythresh)
  cpp1 <- cppicka(xi, obswin, setcov_boundarythresh = setcov_boundarythresh)
  phat <- coverageprob(xi, obswin)
  
  pclns <- paircorr.cvchat(cvchat, cpp1, phat, modifications = modifications, drop = drop) 
  
  return(pclns)
}

#' @describeIn paircorr Generates pair-correlation estimates from
#'   a traditional estimate of covariance, Picka's reduced window estimate of coverage probability,
#'   and the traditional estimate of coverage probability.
#'   If these estimates already exist then \code{paircorr.cvchat} can save significant computation time.
paircorr.cvchat <- function(cvchat, cpp1 = NULL, phat = NULL, modifications = "all", drop = FALSE){
  harmonised <- harmonise.im(cvchat = cvchat, cpp1 = cpp1)
  cvchat <- harmonised$cvchat
  cpp1 <- harmonised$cpp1
  fcns <- list(
         none = pcln_none,
         symm = pcln_symm,
         mattfeldt = pcln_mattfeldt,
         pickaint = pcln_picka_intr,
         pickaH = pcln_picka_H
  )
  if ((modifications == "all")[[1]]) {modifications <- names(fcns)}
  fcnstouse <- fcns[names(fcns) %in% modifications]
  isfunction <- unlist(lapply(modifications, function(x) "function" %in% class(x)))
  modificationsnotused <- modifications[!( (modifications %in% names(fcns)) | isfunction)]
  
  fcnstouse <- c(fcnstouse, modifications[isfunction]) #add user specified modification
  
  if(length(modificationsnotused) > 0){stop(
    paste("The following modifications are not recognised as existing function names or as a function:", modificationsnotused))}
  pclns <- lapply(fcnstouse, function(x) do.call(x, args = list(cvchat = cvchat, cpp1 = cpp1, phat = phat)))
  if (drop & (length(pclns) == 1)) {return(pclns[[1]])}
  else { return(as.imlist(pclns)) }
}


pcln_none <- function(cvchat, cpp1 = NULL, phat = NULL){
  return(cvchat / (phat^2)) 
}

pcln_symm <- function(cvchat, cpp1 = NULL, phat = NULL){
  return((cvchat + reflect.im(cvchat)) / (2 * phat^2)) 
}

pcln_mattfeldt <- function(cvchat, cpp1, phat = NULL){
  return(4 * cvchat /((cpp1 + reflect.im(cpp1))^2) )  
}

pcln_picka_intr <- function(cvchat, cpp1, phat = NULL){
  return(cvchat / (cpp1*reflect.im(cpp1))) 
}

pcln_picka_H <- function(cvchat, cpp1, phat){
  hajek <- phat * (cpp1 + reflect.im(cpp1) - 2 * phat )
  return((cvchat - hajek) / (phat^2))
}

