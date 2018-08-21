#' @title Centred covariance estimation
#' @export cencovariance  cencovariance.cvchat
#' @description 
#' This function estimates the centred covariance of a stationary RACS. 
#' The traditional centred covariance estimator, two 'balanced' estimators suggested by Picka
#'  and a third 'balanced' estimator inspired by one of Picka's pair-correlation estimators.
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


#' @return If \code{drop = TRUE} and a single estimator requested then a
#'   \code{im} object containing the centred covariance estimate. Otherwise a
#'   named \code{imlist} of \code{im} objects containing the centred covariance
#'   estimates for each requested estimator.
#'
#' @keywords spatial nonparametric
#' @details The centred covariance of a stationary RACS is \deqn{\kappa(v) =
#'   C(v) - p^2.}
#'
#'   The estimators available are (see [Chapter 4, hingee2019thesis] for
#'   information): 
#'   \itemize{ 
#'   \item{\code{trad}} the traditional centred
#'   covariance estimator 
#'   \item{\code{mattfeldt}} an estimator inspired by an
#'   `intrinsically' balanced pair-correlation estimator from Picka that was
#'   later studied in an isotropic situation by Mattfeldt and Stoyan
#'   [mattfeldt2000im] 
#'   \item{\code{pickaint}} Picka's intrinsically' balanced
#'   centred covariance estimator [picka2000va]. 
#'   \item{\code{pickaH}} Picka's
#'   additively' balanced centred covariance estimator [picka2000va].
#'   }
#'
#'   Currently computes centred covariance using \code{\link{racscovariance}}.
#'
#' @examples
#' xi <- heather$coarse
#' obswin <- Frame(xi)
#' cencovariance(xi, obswin, modifications = "all")
#' 
#' @describeIn cencovariance Centred covariance estimates from a binary map.
cencovariance <- function(xi, obswin = NULL,
        setcov_boundarythresh = NULL,
        modifications = "all",
        drop = FALSE){
  cvchat <- tradcovarest(xi, obswin, setcov_boundarythresh = setcov_boundarythresh)
  cpp1 <- cppicka(xi, obswin, setcov_boundarythresh = setcov_boundarythresh)
  phat <- coverageprob(xi, obswin)
  
  ccvchats <- cencovariance.cvchat(cvchat, cpp1, phat, modifications = modifications, drop = drop) 
  
  return(ccvchats)
}

#' @describeIn cencovariance Generates centred covariances estimates from
#'   a traditional estimate of covariance, Picka's reduced window estimate of coverage probability,
#'   and the traditional estimate of coverage probability.
#'   If these estimates already exist then \code{cencovariance.cvchat} can save significant computation time.
cencovariance.cvchat <- function(cvchat, cpp1 = NULL, phat = NULL,
        setcov_boundarythresh = NULL,
        modifications = "all",
        drop = FALSE){
  
  cvchats <- racscovariance.cvchat(cvchat, cpp1, phat, modifications = modifications, drop = drop) 
  
  return(cvchats - phat^2)
}
