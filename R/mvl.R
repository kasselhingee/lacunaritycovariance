#' @title Mass variance lacunarity estimatation using all estimators
#' @export mvl mvl.cvchat
#' @description Estimates mass variance lacunarity (MVL) using all estimators described in [Hingee2019**] from binary maps for square box.
#' It calls the functions \code{mvlc}, \code{mvlg}, \code{mvlcc} and \code{mvlgb}.

#' @param xiim A \pkg{spatstat} \code{im} object with pixel values that are either TRUE, FALSE or NA. TRUE represents foreground, FALSE respresents background and NA represents unobserved locations.
#' @param boxwidths A list of box boxwidths
#' @param estimators A list of estimator names - see details for possibilities. \code{estimators = "all"} will select all estimators.
#' @param includenormed A logical value. If TRUE then MVL estimates normalised by the MVL values at zero will be included in a returned list of fv objects
#' @param setcov_boundarythresh Any vector \eqn{v} such that set covariance of the observation window is smaller than this threshold
#' is given a covariance estimate (and other similar estimate) of NA to avoid instabilities caused by dividing by very small areas.
#' If NULL is supplied (default) then 1E-6 is used.
#' @param phat Traditional estimate of coverage probability.
#' @param cvchat Traditional estimate of covariance (often from \code{tradcovarest}).
#' @param cpp1 Picka's estimate of coverage probability (often from \code{cppicka}).

#' @return An \code{fv} object.

#' @details
#' The function is not able to estimate MVL for non-square boxes as the gliding box estimator is included.
#' To estimate MVL for non-square boxes use \code{mvlcc} or \code{mvlg} directly.
#' 
#' The estimators available are
#' \itemize{
#' \item{\code{"MVLc"}} The unmodified (unbalanced) covariance estimator provided by \code{\link{mvlc}}
#' \item{\code{"MVLgb"}} The Gliding-Box estimator of Alain and Cloitre [allain1991ch]. Calls \code{\link{mvlgb}}
#' \item{\code{"MVLg.mattfeldt"}} See help for \code{\link{mvlg}}
#' \item{\code{"MVLg.pickaint"}} See help for \code{\link{mvlg}}
#' \item{\code{"MVLg.pickaH"}} See help for \code{\link{mvlg}}
#' \item{\code{"MVLcc.mattfeldt"}} See help for \code{\link{mvlcc}}
#' \item{\code{"MVLcc.pickaint"}} See help for \code{\link{mvlcc}}
#' \item{\code{"MVLcc.pickaH"}} See help for \code{\link{mvlcc}}
#' }

#' @examples 
#' xi <- heather$coarse
#' xiim <- as.im(xi, value = TRUE, na.replace = FALSE)
#' mvlests <- mvl(xiim, seq(1, 10, by = 0.1))

#' @describeIn mvl Computes MVL estimates from a binary map.
mvl <- function(xiim, boxwidths,
                           estimators = c("MVLg.mattfeldt", "MVLg.pickaint", "MVLg.pickaH",
                                          "MVLcc.mattfeldt", "MVLcc.pickaint",
                                          "MVLc", "MVLgb"),
                includenormed = FALSE,
                setcov_boundarythresh = 1E-6){
  if ("all" %in% estimators){
    estimators = c("MVLg.mattfeldt", "MVLg.pickaint", "MVLg.pickaH",
     "MVLcc.mattfeldt", "MVLcc.pickaint",
     "MVLc", "MVLgb")
  }
  mvlgestimaterequests <- estimators %in% MVLgestimatornames
  mvlccestimaterequests <- estimators %in% MVLccestimatornames
  
  mvl.ests <- list()
  cpp1 <- cvchat <- NULL
  mvlcovarbased <- mvlgb.est <- NULL
  
  phat <- coverageprob(xiim)
  if(sum(mvlgestimaterequests) + sum(mvlccestimaterequests) + ("MVLc" %in% estimators) > 0){
    cvchat <- tradcovarest(xiim, setcov_boundarythresh = setcov_boundarythresh)
  }
  if (sum(mvlgestimaterequests) + sum(mvlccestimaterequests) > 0){
    cpp1 <- cppicka(xiim, setcov_boundarythresh = setcov_boundarythresh)
  }
  if(sum(mvlgestimaterequests) + sum(mvlccestimaterequests) + ("MVLc" %in% estimators) > 0){
    #function that computes the covariance-based estimates of MVL
    mvlcovarbased <- mvl.cvchat(boxwidths = boxwidths, estimators = estimators, phat = phat, cvchat = cvchat, cpp1 = cpp1)
    mvl.ests <- c(mvl.ests, mvlcovarbased)
  }
  
  #the MVLgb estimate
  if ("MVLgb" %in% estimators){
    mvlgb.est <- mvlgb(sidelengths = boxwidths, xiim = xiim)
    if (sum(!vapply(mvlgb.est[,fvnames(mvlgb.est), drop = TRUE], is.na, FUN.VALUE = TRUE)) < 2){
      warning("mvlgb() returns estimates for 1 or fewer of the provided box widths. Results from mvlgb() will be ignored from the final results.")
      mvlgb.est <- NULL
    }
    mvl.ests <- c(mvl.ests, list(mvlgb = mvlgb.est))
  }
  
  mvl.ests <- mvl.ests[!vapply(mvl.ests, is.null, FUN.VALUE = FALSE)]
  if (any(!vapply(mvl.ests[-1], function(x) compatible.fv(A = mvl.ests[[1]], B = x), FUN.VALUE = FALSE))){
    warning("Some MVL estimates have differing argument values. These will be harmonised.")
    mvl.ests <- harmonise.fv(mvl.ests)
  }
  mvls.fv <- collapse.fv(mvl.ests, different = "MVL")
  names(mvls.fv) <- c(fvnames(mvls.fv, ".x"), names(mvl.ests))
  
  allfvs <- list(mvl.est = mvls.fv)
  
  if (includenormed){
    #compute MVLs normalised at zero
    normdmvls <- eval.fv(mvls.fv / ( phat * (1 - phat) / phat^2), relabel = FALSE)
    normdmvls <- prefixfv(normdmvls,
                     tagprefix="n_",
                     descprefix="normalised at zero",
                     lablprefix="plain(nrmd)~")
    allfvs <- c(allfvs, list(normdmvls = normdmvls))
  }
  
  return(allfvs)
}

#' @describeIn mvl Computes covariance-based estimator of MVL from the traditional estimate of covariance,
#'  Picka's reduced window coverage probability estimates and the traditional coverage probability estimate.
mvl.cvchat <- function(boxwidths,
                      estimators = c("MVLg.mattfeldt", "MVLg.pickaint", "MVLg.pickaH",
                                     "MVLcc.mattfeldt", "MVLcc.pickaint",
                                     "MVLc"),
                      phat = NULL,
                      cvchat = NULL,
                      cpp1 = NULL){
  if ("all" %in% estimators){
    estimators = c("MVLg.mattfeldt", "MVLg.pickaint", "MVLg.pickaH",
     "MVLcc.mattfeldt", "MVLcc.pickaint",
     "MVLc")
  }
  mvlgestimaterequests <- estimators %in% MVLgestimatornames
  mvlccestimaterequests <- estimators %in% MVLccestimatornames
  mvlgs <- mvlccs <- mvlc.est <- mvlgb.est <- NULL
  
  if (sum(mvlgestimaterequests) + sum(mvlccestimaterequests) > 0){
    stopifnot(!is.null(cpp1), !is.null(cvchat), !is.null(phat))
    if (sum(mvlgestimaterequests) > 0){
      pcln.ests <- paircorr.cvchat(cvchat, cpp1 = cpp1, phat = phat, modifications = gsub("MVLg.", "", estimators[mvlgestimaterequests]), drop = FALSE)
      mvlgs <- lapply(pcln.ests, FUN = mvlg, boxes = boxwidths)
    }
    if (sum(mvlccestimaterequests) > 0){
      ccvc.ests <- cencovariance.cvchat(cvchat, cpp1, phat, modifications = gsub("MVLcc.", "", estimators[mvlccestimaterequests]), drop = FALSE)
      mvlccs <- lapply(ccvc.ests, FUN = mvlcc, p = phat, boxes = boxwidths)
    }
  }
  if ("MVLc" %in% estimators){
    stopifnot(!is.null(cvchat), !is.null(phat))
    mvlc.est <- mvlc(boxes = boxwidths, covariance = cvchat, p = phat)
  }
  mvl.ests <- c(mvlg = mvlgs, mvlcc = mvlccs, list(mvlc = mvlc.est))
  mvl.ests <- mvl.ests[!vapply(mvl.ests, is.null, FUN.VALUE = FALSE)]
  return(mvl.ests)
}



MVLgestimatornames <- c("MVLg.none", "MVLg.mattfeldt", "MVLg.pickaint", "MVLg.pickaH")
MVLccestimatornames <- c("MVLcc.none", "MVLcc.mattfeldt", "MVLcc.pickaint", "MVLcc.pickaH")
