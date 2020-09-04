#' @title Centred covariance estimation
#' @export cencovariance  cencovariance.cvchat
#' @description 
#' This function estimates the centred covariance of a stationary RACS. 
#' Available estimators are the plug-in moment centred covariance estimator, two 'balanced' estimators suggested by Picka (2000),
#'  and a third 'balanced' estimator inspired by one of Picka's pair-correlation estimators.
#' @author{Kassel Liam Hingee}

#' @param xi An observation of a RACS of interest as a full binary map (as an \code{im} object) or as the foreground set (as an \code{owin} object).
#' In the latter case the observation window, \code{obswin}, must be supplied.
#' @param obswin If \code{xi} is an \code{owin} object then \code{obswin} is an
#'   \code{owin} object that specifies the observation window.
#' @param setcov_boundarythresh To avoid instabilities caused by dividing by very small quantities, if the set covariance of the observation window
#'  is smaller than \code{setcov_boundarythresh}, then the covariance is given a value of NA. 
#' @param phat The usual estimate of coverage probability,
#'  which is the observed foreground area in \code{xi} divided by the total area of the observation window.
#'  See \code{\link{coverageprob}} for more information.
#' @param cvchat The plug-in moment estimate of covariance as an \code{im} object. 
#' Typically created with \code{\link{plugincvc}}.
#' @param cpp1 Picka's reduced window estimate of coverage probability as an \code{im} object - used in improved (balanced) covariance estimators.
#' Can be generated using \code{\link{cppicka}}.
#' @param estimators A list of strings specifying estimators to use. 
#' See details.
#' \code{estimators = "all"} will select all available estimators.  
#' @param drop If TRUE and one estimator selected then the returned value will be a single \code{im} object and not a list of \code{im} object.


#' @return If \code{drop = TRUE} and only one estimator is requested then a
#'   \code{im} object containing the centred covariance estimate is returned. Otherwise a
#'   named \code{imlist} of \code{im} objects containing the centred covariance
#'   estimates for each requested estimator.
#'
#' @keywords spatial nonparametric
#' @details The centred covariance of a stationary RACS is \deqn{\kappa(v) =
#'   C(v) - p^2.}
#'
#'   The estimators available are (see (Section 3.4, Hingee, 2019) for
#'   more information): 
#'   \itemize{ 
#'   \item{\code{plugin}} the plug-in moment centred
#'   covariance estimator 
#'   \item{\code{mattfeldt}} an estimator inspired by an
#'   'intrinsically' balanced pair-correlation estimator from Picka (1997) that was
#'   later studied in an isotropic situation by Mattfeldt and Stoyan
#'   (Mattfeldt and Stoyan, 2000) 
#'   \item{\code{pickaint}} Picka's 'intrinsically' balanced
#'   centred covariance estimator (Picka, 2000). 
#'   \item{\code{pickaH}} Picka's
#'   'additively' balanced centred covariance estimator (Picka, 2000).
#'   }
#'
#'   Currently computes centred covariance using \code{\link{racscovariance}}.
#'
#' @references
#' Hingee, K.L. (2019) \emph{Spatial Statistics of Random Closed Sets for Earth Observations}. PhD: Perth, Western Australia: University of Western Australia. Submitted.
#'
#' Mattfeldt, T. and Stoyan, D. (2000) Improved estimation of the pair correlation function of random sets. \emph{Journal of Microscopy}, 200, 158-173.
#'
#' Picka, J.D. (1997) \emph{Variance-Reducing Modifications for Estimators of Dependence in Random Sets}. Ph.D.: Illinois, USA: The University of Chicago.
#' 
#' Picka, J.D. (2000) Variance reducing modifications for estimators of standardized moments of random sets. \emph{Advances in Applied Probability}, 32, 682-700.
#'
#' @examples
#' xi <- heather$coarse
#' obswin <- Frame(xi)
#' cencovariance(xi, obswin, estimators = "all")
#' 
#' @describeIn cencovariance Centred covariance estimates from a binary map.
cencovariance <- function(xi, obswin = NULL,
        setcov_boundarythresh = NULL,
        estimators = "all",
        drop = FALSE){
  cvchat <- plugincvc(xi, obswin, setcov_boundarythresh = setcov_boundarythresh)
  cpp1 <- cppicka(xi, obswin, setcov_boundarythresh = setcov_boundarythresh)
  phat <- coverageprob(xi, obswin)
  
  ccvchats <- cencovariance.cvchat(cvchat, cpp1, phat, estimators = estimators, drop = drop) 
  
  return(ccvchats)
}

#' @describeIn cencovariance Generates centred covariances estimates from
#'   a plug-in moment estimate of covariance, Picka's reduced window estimate of coverage probability,
#'   and the plug-in moment estimate of coverage probability.
#'   If these estimates already exist, then \code{\link{cencovariance.cvchat}} saves significant computation time over \code{cencovariance}.
cencovariance.cvchat <- function(cvchat, cpp1 = NULL, phat = NULL,
        setcov_boundarythresh = NULL,
        estimators = "all",
        drop = FALSE){
  
  cvchats <- racscovariance.cvchat(cvchat, cpp1, phat, estimators = estimators, drop = drop) 
  
  return(cvchats - phat^2)
}
