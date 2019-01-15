#' @title Gliding box lacunarity estimatation using all estimators
#' @export gbl gbl.cvchat
#' @description Estimates gliding box lacunarity (GBL) using all estimators described in (Hingee et al., 2017) from binary maps for square boxes.
#' It calls the functions \code{gblc}, \code{gblg}, \code{gblcc} and \code{gbltrad}.

#' @param xi An observation of a RACS of interest as a full binary map (in \code{im} format) or as the foreground set (in \code{owin} format).
#' In the latter case the observation window, \code{obswin}, must be supplied.
#' See \code{\link{stationaryracsinference-package}} for details.
#' @param obswin If \code{xi} is an \code{owin} object then \code{obswin} is an
#'   \code{owin} object that specifies the observation window.
#' @param boxwidths A list of box widths
#' @param estimators A list of estimator names - see details for possibilities. \code{estimators = "all"} will select all estimators.
#' @param includenormed A logical value. If TRUE then GBL estimates normalised by the GBL values at zero will be included in a returned list of fv objects
#' @param setcov_boundarythresh Any vector \eqn{v} such that set covariance of the observation window is smaller than this threshold
#' is given a covariance estimate (and other similar estimate) of NA to avoid instabilities caused by dividing by very small areas.
#' If NULL is supplied (default) then 1E-6 is used.
#' @param phat Traditional estimate of coverage probability.
#' @param cvchat Traditional estimate of covariance (often from \code{tradcovarest}).
#' @param cpp1 Picka's estimate of coverage probability (often from \code{cppicka}).

#' @return An \code{fv} object.

#' @details
#' The function is not able to estimate GBL for non-square boxes as the gliding box estimator is included.
#' To estimate GBL for non-square boxes use \code{gblcc} or \code{gblg} directly.
#' 
#' If \code{xi} is in \code{owin} format then \code{obswin} and \code{xi} are converted
#'  into a binary map in \code{im} format using \code{\link[spatstat]{as.im}}
#' 
#' The estimators available are
#' \itemize{
#' \item{\code{"GBLc"}} The unmodified (unbalanced) covariance estimator provided by \code{\link{gblc}}
#' \item{\code{"GBLgb"}} The Gliding-Box estimator of Allain and Cloitre (1991). Calls \code{\link{gbltrad}}
#' \item{\code{"GBLg.mattfeldt"}} See help for \code{\link{gblg}}
#' \item{\code{"GBLg.pickaint"}} See help for \code{\link{gblg}}
#' \item{\code{"GBLg.pickaH"}} See help for \code{\link{gblg}}
#' \item{\code{"GBLcc.mattfeldt"}} See help for \code{\link{gblcc}}
#' \item{\code{"GBLcc.pickaint"}} See help for \code{\link{gblcc}}
#' \item{\code{"GBLcc.pickaH"}} See help for \code{\link{gblcc}}
#' }

#' @references
#' Allain, C. and Cloitre, M. (1991) Characterizing the lacunarity of random and deterministic fractal sets. \emph{Physical Review A}, 44, 3552-3558.
#' 
#' Hingee K, Baddeley A, Caccetta P, Nair G (2017). Computation of lacunarity from covariance of spatial binary maps. \emph{Journal of Agricultural, Biological and Environmental Statistics}. Submitted.


#' @examples 
#' xi <- heather$coarse
#' xi <- as.im(xi, value = TRUE, na.replace = FALSE)
#' gblests <- gbl(xi, seq(1, 10, by = 0.1))

#' @describeIn gbl Computes GBL estimates from a binary map.
gbl <- function(xi, boxwidths,
                           estimators = c("GBLg.mattfeldt", "GBLg.pickaint", "GBLg.pickaH",
                                          "GBLcc.mattfeldt", "GBLcc.pickaint",
                                          "GBLc", "GBLgb"),
                obswin = NULL,
                includenormed = FALSE,
                setcov_boundarythresh = 1E-6){
  if ("all" %in% estimators){
    estimators = c("GBLg.mattfeldt", "GBLg.pickaint", "GBLg.pickaH",
     "GBLcc.mattfeldt", "GBLcc.pickaint", "GBLcc.pickaH",
     "GBLc", "GBLgb")
  }
  gblgestimaterequests <- estimators %in% GBLgestimatornames
  gblccestimaterequests <- estimators %in% GBLccestimatornames
  
  gbl.ests <- list()
  cpp1 <- cvchat <- NULL
  gblcovarbased <- gbltrad.est <- NULL
  
  #if an owin is passed convert to an im object 
  if (is.owin(xi)) {
    if (is.null(obswin)) {stop("obswin must be provide if xi is owin")}
    xi <- as.im(xi, value = TRUE, na.replace = FALSE)
    xi[setminus.owin(Frame(xi), obswin)] <- NA
  }
  stopifnot(isbinarymap(xi))
  
  phat <- coverageprob(xi)
  if(sum(gblgestimaterequests) + sum(gblccestimaterequests) + ("GBLc" %in% estimators) > 0){
    cvchat <- tradcovarest(xi, setcov_boundarythresh = setcov_boundarythresh)
  }
  if (sum(gblgestimaterequests) + sum(gblccestimaterequests) > 0){
    cpp1 <- cppicka(xi, setcov_boundarythresh = setcov_boundarythresh)
  }
  if(sum(gblgestimaterequests) + sum(gblccestimaterequests) + ("GBLc" %in% estimators) > 0){
    #function that computes the covariance-based estimates of GBL
    gblcovarbased <- gbl.cvchat(boxwidths = boxwidths, estimators = estimators, phat = phat, cvchat = cvchat, cpp1 = cpp1)
    gbl.ests <- c(gbl.ests, gblcovarbased)
  }
  
  #the GBLgb estimate
  if ("GBLgb" %in% estimators){
    gbltrad.est <- gbltrad(boxwidths = boxwidths, xiim = xi)
    if (sum(!vapply(gbltrad.est[,fvnames(gbltrad.est), drop = TRUE], is.na, FUN.VALUE = TRUE)) < 2){
      warning("gbltrad() returns estimates for 1 or fewer of the provided box widths. Results from gbltrad() will be ignored from the final results.")
      gbltrad.est <- NULL
    }
    gbl.ests <- c(gbl.ests, list(gbltrad = gbltrad.est))
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

#' @describeIn gbl Computes covariance-based estimator of GBL from the traditional estimate of covariance,
#'  Picka's reduced window coverage probability estimates and the traditional coverage probability estimate.
gbl.cvchat <- function(boxwidths,
                      estimators = c("GBLg.mattfeldt", "GBLg.pickaint", "GBLg.pickaH",
                                     "GBLcc.mattfeldt", "GBLcc.pickaint",
                                     "GBLc"),
                      phat = NULL,
                      cvchat = NULL,
                      cpp1 = NULL){
  if ("all" %in% estimators){
    estimators = c("GBLg.mattfeldt", "GBLg.pickaint", "GBLg.pickaH",
     "GBLcc.mattfeldt", "GBLcc.pickaint",
     "GBLc")
  }
  gblgestimaterequests <- estimators %in% GBLgestimatornames
  gblccestimaterequests <- estimators %in% GBLccestimatornames
  gblgs <- gblccs <- gblc.est <- gbltrad.est <- NULL
  
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



GBLgestimatornames <- c("GBLg.trad", "GBLg.mattfeldt", "GBLg.pickaint", "GBLg.pickaH")
GBLccestimatornames <- c("GBLcc.trad", "GBLcc.mattfeldt", "GBLcc.pickaint", "GBLcc.pickaH")
