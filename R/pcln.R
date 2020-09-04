#' @title Balanced estimation of pair-correlation.
#' @export paircorr  paircorr.cvchat
#' @description 
#' Estimates the pair-correlation function of a stationary RACS. 
#' The plug-in moment pair-correlation estimator and three 'balanced' estimators suggested by Picka (2000)
#' are available.
#' @author{Kassel Liam Hingee}

#' @param xi An observation of a RACS of interest as a full binary map (as an \code{im} object) or as the foreground set (as an \code{owin} object).
#' In the latter case the observation window, \code{obswin}, must be supplied.
#' See \code{\link{lacunaritycovariance-package}} for details.
#' @param obswin If \code{xi} is an \code{owin} object then \code{obswin} is an
#'   \code{owin} object that specifies the observation window.
#' @param setcov_boundarythresh To avoid instabilities caused by dividing by very small quantities, if the set covariance of the observation window
#'  is smaller than \code{setcov_boundarythresh}, then the covariance is given a value of NA. 
#' @param phat The plug-in moment estimate of coverage probability,
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
#'  \code{estimators = "all"} will select all inbuilt estimators. See details. 

#' @return If \code{drop = TRUE} and a single estimator is requested then an
#'   \code{im} object containing the pair-correlation estimate is returned. Otherwise a
#'   named \code{imlist} of \code{im} objects containing the pair-correlation
#'   estimates for each requested estimator.
#'
#'   


#' @keywords spatial nonparametric
#' @details The pair-correlation of a stationary RACS is 
#' \deqn{g(v) = C(v) / p^2.}
#'
#'   The estimators available are (see (Hingee, 2019) for
#'   more information): 
#'   \itemize{ 
#'   \item{\code{plugin}} the plug-in moment pair-correlation estimator which is \eqn{Chat(v) / (phat^2)}, where \eqn{Chat} and \eqn{phat} are 
#' the plug-in moment estimate of covariance and the usual estimate of coverage probability, respectively.
#'   \item{\code{mattfeldt}} an 'intrinsically' balanced pair-correlation estimator suggested by Picka (1997).
#'   A similar isotropic pair-correlation estimator was later studied by Mattfeldt and Stoyan (2000).
#'   \item{\code{pickaint}} Picka's 'intrinsically' balanced pair-correlation estimator (Picka, 2000). 
#'   \item{\code{pickaH}} Picka's 'additively' balanced pair-correlation estimator (Picka, 2000).
#'   }
#'
#' @references
#' Hingee, K.L. (2019) \emph{Spatial Statistics of Random Closed Sets for Earth Observations}. PhD: Perth, Western Australia: University of Western Australia. Submitted.
#' 
#' Mattfeldt, T. and Stoyan, D. (2000) Improved estimation of the pair correlation function of random sets. \emph{Journal of Microscopy}, 200, 158-173.
#' 
#' Picka, J.D. (1997) \emph{Variance-Reducing Modifications for Estimators of Dependence in Random Sets}. Ph.D.: Illinois, USA: The University of Chicago.
#' 
#' Picka, J.D. (2000) Variance reducing modifications for estimators of standardized moments of random sets. \emph{Advances in Applied Probability}, 32, 682-700.



#' @examples
#' xi <- as.im(heather$coarse, na.replace = 0, eps = 4 * heather$coarse$xstep)
#'
#' # Estimate pair correlation from a binary map
#' pclns_directest <- paircorr(xi, estimators = "all")
#' 
#' phat <- coverageprob(xi)
#' cvchat <- plugincvc(xi)
#' cpp1 <- cppicka(xi)
#' 
#' # Compute pair correlation estimates from estimates covariance,
#' # coverage probability and Picka's reduced-window coverage probability.
#' pclns_fromcvc <- paircorr.cvchat(cvchat, cpp1, phat, estimators = "all")

#' @describeIn paircorr Estimates pair-correlation from a binary map.
paircorr <- function(xi, obswin = NULL,
                  setcov_boundarythresh = NULL,
                  estimators = "all",
                  drop = FALSE){
  cvchat <- plugincvc(xi, obswin, setcov_boundarythresh = setcov_boundarythresh)
  cpp1 <- cppicka(xi, obswin, setcov_boundarythresh = setcov_boundarythresh)
  phat <- coverageprob(xi, obswin)
  
  pclns <- paircorr.cvchat(cvchat, cpp1, phat, estimators = estimators, drop = drop) 
  
  return(pclns)
}

#' @describeIn paircorr Generates pair-correlation estimates from
#'   the plug-in moment estimates of covariance, Picka's reduced window estimate of coverage probability,
#'   and the coverage fraction (which is an unbiased estimate of the coverage probability).
#'   If these estimates already exist then \code{paircorr.cvchat} can save significant computation time.
paircorr.cvchat <- function(cvchat, cpp1 = NULL, phat = NULL, estimators = "all", drop = FALSE){
  harmonised <- harmonise.im(cvchat = cvchat, cpp1 = cpp1)
  cvchat <- harmonised$cvchat
  cpp1 <- harmonised$cpp1
  fcns <- list(
         plugin = pcln_plugin,
         symm = pcln_symm,
         mattfeldt = pcln_mattfeldt,
         pickaint = pcln_picka_intr,
         pickaH = pcln_picka_H
  )
  if ((estimators == "all")[[1]]) {estimators <- names(fcns)}
  fcnstouse <- fcns[names(fcns) %in% estimators]
  isfunction <- unlist(lapply(estimators, function(x) "function" %in% class(x)))
  estimatorsnotused <- estimators[!( (estimators %in% names(fcns)) | isfunction)]
  
  fcnstouse <- c(fcnstouse, estimators[isfunction]) #add user specified estimator
  
  if(length(estimatorsnotused) > 0){stop(
    paste("The following estimators are not recognised as existing function names or as a function:", estimatorsnotused))}
  pclns <- lapply(fcnstouse, function(x) do.call(x, args = list(cvchat = cvchat, cpp1 = cpp1, phat = phat)))
  if (drop & (length(pclns) == 1)) {return(pclns[[1]])}
  else { return(as.imlist(pclns)) }
}


pcln_plugin <- function(cvchat, cpp1 = NULL, phat = NULL){
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

