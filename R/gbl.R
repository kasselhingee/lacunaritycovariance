#' @title Gliding box lacunarity estimation using all estimators
#' @export gbl gbl.cvchat
#' @description Estimates gliding box lacunarity (GBL) for square boxes using all estimators described in (Hingee et al., 2017).
#' It calls the functions \code{\link{gblc}}, \code{\link{gblg}}, \code{\link{gblcc}} and \code{\link{gblemp}}.

#' @param xi An observation of a RACS of interest as a full binary map (as an \code{im} object) or as the foreground set (as an \code{owin} object).
#' In the latter case the observation window, \code{obswin}, must be supplied.
#' @param obswin If \code{xi} is an \code{owin} object then \code{obswin} is an
#'   \code{owin} object that specifies the observation window.
#' @param boxwidths A list of box widths
#' @param estimators A vector of estimator names - see details for possible names. \code{estimators = "all"} will select all estimators.
#' @param includenormed A logical value. If TRUE then GBL estimates normalised by the GBL values at zero will be included in a returned list of \code{fv} objects
#' @param setcov_boundarythresh To avoid instabilities caused by dividing by very small quantities, if the set covariance of the observation window
#'  is smaller than \code{setcov_boundarythresh}, then the covariance is given a value of NA. 
#' If NULL is supplied (default) then \code{setcov_boundarythresh} is set to 1E-6.
#' @param phat  The fraction of foreground area in the observation window, which is the usual estimator of coverage probability given by \code{\link{coverageprob}}.
#' @param cvchat The plug-in moment covariance estimate (often from \code{\link{plugincvc}}).
#' @param cpp1 Picka's estimate of coverage probability for subsets of the observation window. See \code{\link{cppicka}}.

#' @return An \code{fv} object.

#' @details
#' As empirical GBL is one of the GBL estimators available through this function, non-square boxes are not allowed.
#' To estimate GBL for non-square boxes use \code{gblcc} or \code{gblg} directly.
#' 
#' If \code{xi} is an \code{owin} object then \code{obswin} and \code{xi} are converted
#'  into an \code{im} object using \code{\link[spatstat]{as.im}}
#' 
#' The estimators available are
#' \itemize{
#' \item{\code{"GBLc"}} The unmodified (unbalanced) covariance estimator provided by \code{\link{gblc}}
#' \item{\code{"GBLemp"}} Empirical gliding box lacunarity (Allain and Cloitre, 1991). Calls \code{\link{gblemp}}
#' \item{\code{"GBLg.mattfeldt"}} See help for \code{\link{gblg}} and \code{\link{paircorr}}
#' \item{\code{"GBLg.pickaint"}} See help for \code{\link{gblg}} and \code{\link{paircorr}}
#' \item{\code{"GBLg.pickaH"}} See help for \code{\link{gblg}} and \code{\link{paircorr}}
#' \item{\code{"GBLcc.mattfeldt"}} See help for \code{\link{gblcc}}
#' \item{\code{"GBLcc.pickaint"}} See help for \code{\link{gblcc}}
#' \item{\code{"GBLcc.pickaH"}} See help for \code{\link{gblcc}}
#' }
#' 
#' The default, GBLcc.pickaH, is a method based on centred covariance.
#' Centred covariance approaches zero for large vectors, and are thus easier to integrate with the set covariance of the boxes.
#'
#' The set covariance of the boxes is computed empirically using \pkg{spatstat}'s \code{\link[spatstat]{setcov}} function, which converts \eqn{B} into a binary pixel mask using \code{\link[spatstat]{as.mask}} defaults. Computation speed can be increased by setting a small default number of pixels, \code{npixel}, in \pkg{spatstat}'s global options (accessed through \code{\link[spatstat]{spatstat.options}}), however fewer pixels also decreases the accuracy of the GBL computation.

#' @references
#' Allain, C. and Cloitre, M. (1991) Characterizing the lacunarity of random and deterministic fractal sets. \emph{Physical Review A}, 44, 3552-3558.
#' 
#' Hingee K, Baddeley A, Caccetta P, Nair G (2019). Computation of lacunarity from covariance of spatial binary maps. \emph{Journal of Agricultural, Biological and Environmental Statistics}, 24, 264-288. DOI: 10.1007/s13253-019-00351-9.


#' @examples 
#' xi <- as.im(heather$coarse, value = TRUE,
#'             na.replace = FALSE)
#' boxwidths <- seq(1, 10, by = 0.5)
#' 
#' # reduce resolution in setcov() for faster (less accurate) computation 
#' oldopt <- spatstat.options()
#' spatstat.options("npixel" = 2^5)
#' 
#' gblests <- gbl(xi, boxwidths, estimators = "GBLcc.pickaH")
#' spatstat.options(oldopt)


#' @describeIn gbl computes GBL estimates from a binary map.
gbl <- function(xi, boxwidths,
                           estimators = c("GBLcc.pickaH"),
                obswin = NULL,
                includenormed = FALSE,
                setcov_boundarythresh = 1E-6){
  if ("all" %in% estimators){
    estimators = c(GBLgestimatornames, GBLccestimatornames, 
     "GBLc", "GBLemp")
  }
  gblgestimaterequests <- estimators %in% GBLgestimatornames
  gblccestimaterequests <- estimators %in% GBLccestimatornames
  
  gbl.ests <- list()
  cpp1 <- cvchat <- NULL
  gblcovarbased <- gblemp.est <- NULL
  
  #if an owin is passed convert to an im object 
  if (is.owin(xi)) {
    if (is.null(obswin)) {stop("obswin must be provide if xi is owin")}
    xi <- as.im(xi, value = TRUE, na.replace = FALSE)
    xi[setminus.owin(Frame(xi), obswin)] <- NA
  }
  stopifnot(isbinarymap(xi))
  
  phat <- coverageprob(xi)
  if(sum(gblgestimaterequests) + sum(gblccestimaterequests) + ("GBLc" %in% estimators) > 0){
    cvchat <- plugincvc(xi, setcov_boundarythresh = setcov_boundarythresh)
  }
  if (sum(gblgestimaterequests) + sum(gblccestimaterequests) > 0){
    cpp1 <- cppicka(xi, setcov_boundarythresh = setcov_boundarythresh)
  }
  if(sum(gblgestimaterequests) + sum(gblccestimaterequests) + ("GBLc" %in% estimators) > 0){
    #function that computes the covariance-based estimates of GBL
    gblcovarbased <- gbl.cvchat(boxwidths = boxwidths, estimators = estimators, phat = phat, cvchat = cvchat, cpp1 = cpp1)
    gbl.ests <- c(gbl.ests, gblcovarbased)
  }
  
  #the GBLemp estimate
  if ("GBLemp" %in% estimators){
    gblemp.est <- gblemp(boxwidths = boxwidths, xiim = xi)
    if (sum(!vapply(gblemp.est[,fvnames(gblemp.est), drop = TRUE], is.na, FUN.VALUE = TRUE)) < 2){
      warning("gblemp() returns estimates for 1 or fewer of the provided box widths. Results from gblemp() will be ignored from the final results.")
      gblemp.est <- NULL
    }
    gbl.ests <- c(gbl.ests, list(gblemp = gblemp.est))
  }
  
  gbl.ests <- gbl.ests[!vapply(gbl.ests, is.null, FUN.VALUE = FALSE)]
  if (any(!vapply(gbl.ests[-1], function(x) compatible.fv(A = gbl.ests[[1]], B = x), FUN.VALUE = FALSE))){
    warning("Some GBL estimates have differing argument values. These will be harmonised.")
    gbl.ests <- harmonise.fv(gbl.ests)
  }
  gbls.fv <- collapse.fv(gbl.ests, different = "GBL")
  names(gbls.fv) <- c(fvnames(gbls.fv, ".x"), names(gbl.ests))
  
  allfvs <- list(gbl.est = gbls.fv)
  
  if (includenormed){
    #compute GBLs normalised at zero
    normdgbls <- eval.fv((gbls.fv - 1)/ ( phat * (1 - phat) / phat^2), relabel = FALSE)
    normdgbls <- prefixfv(normdgbls,
                     tagprefix="n_",
                     descprefix="standardised ",
                     lablprefix="plain(nrmd)~")
    allfvs <- c(allfvs, list(normdgbls = normdgbls))
  }
  
  return(allfvs)
}

#' @describeIn gbl computes covariance-based estimator of GBL from the plug-in moment estimate of covariance,
#'  Picka's reduced window coverage probability estimates (see \code{\link{cppicka}}) and the usual coverage probability estimate, \code{phat}.
gbl.cvchat <- function(boxwidths,
                      estimators = c("GBLg.mattfeldt", "GBLg.pickaint", "GBLg.pickaH",
                                     "GBLcc.mattfeldt", "GBLcc.pickaint",
                                     "GBLc"),
                      phat = NULL,
                      cvchat = NULL,
                      cpp1 = NULL){
  if ("all" %in% estimators){
    estimators = c(GBLgestimatornames, GBLccestimatornames, 
     "GBLc", "GBLemp")
  }
  
  gblgestimaterequests <- estimators %in% GBLgestimatornames
  gblccestimaterequests <- estimators %in% GBLccestimatornames
  gblgs <- gblccs <- gblc.est <- gblemp.est <- NULL
  
  if (sum(gblgestimaterequests) + sum(gblccestimaterequests) > 0){
    stopifnot(!is.null(cpp1), !is.null(cvchat), !is.null(phat))
    if (sum(gblgestimaterequests) > 0){
      pcln.ests <- paircorr.cvchat(cvchat, cpp1 = cpp1, phat = phat, estimators = gsub("GBLg.", "", estimators[gblgestimaterequests]), drop = FALSE)
      gblgs <- lapply(pcln.ests, FUN = gblg, boxes = boxwidths)
    }
    if (sum(gblccestimaterequests) > 0){
      ccvc.ests <- cencovariance.cvchat(cvchat, cpp1, phat, estimators = gsub("GBLcc.", "", estimators[gblccestimaterequests]), drop = FALSE)
      gblccs <- lapply(ccvc.ests, FUN = gblcc, p = phat, boxes = boxwidths)
    }
  }
  if ("GBLc" %in% estimators){
    stopifnot(!is.null(cvchat), !is.null(phat))
    gblc.est <- gblc(boxes = boxwidths, covariance = cvchat, p = phat)
  }
  gbl.ests <- c(gblg = gblgs, gblcc = gblccs, list(gblc = gblc.est))
  gbl.ests <- gbl.ests[!vapply(gbl.ests, is.null, FUN.VALUE = FALSE)]
  return(gbl.ests)
}



GBLgestimatornames <- c("GBLg.plugin", "GBLg.mattfeldt", "GBLg.pickaint", "GBLg.pickaH")
GBLccestimatornames <- c("GBLcc.plugin", "GBLcc.mattfeldt", "GBLcc.pickaint", "GBLcc.pickaH")
