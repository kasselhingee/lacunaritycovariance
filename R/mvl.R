#' @title Compute MVL estimates using all estimators
#' @export mvl
#' @description This function computes all the MVL estimators described in Hingee 2019** for square boxes and raster binary maps.
#' It calls the functions \code{MVLc}, \code{MVLg}, \code{MVLcc} and \code{MVLgb}.

#' @param xiim A \pkg{spatstat} \code{im} object with pixel values that are either TRUE, FALSE or NA. TRUE represents foreground, FALSE respresents background and NA represents unobserved locations.
#' @param boxwidths A list of box boxwidths
#' @param estimators A list of estimator names - see details for possibilities.

#' @return An \code{fv} object.

#' @details
#' The estimators available are
#' \itemize{
#' \item \code{"MVLc"} The unmodified (unbalanced) covariance estimator provided by \code{MVLc()}
#' \item \code{"MVLgb"} The Gliding-Box estimator of Alain and Cloitre **ref. Calls \code{MVLgb()}
#' \item \code{"MVLg.mattfeldt"} See help for \code{MVLg()}
#' \item \code{"MVLg.pickaint"} See help for \code{MVLg()}
#' \item \code{"MVLg.pickaH"} See help for \code{MVLg()}
#' \item \code{"MVLcc.mattfeldt"} See help for \code{MVLcc()}
#' \item \code{"MVLcc.pickaint"} See help for \code{MVLcc()}
#' \item \code{"MVLcc.pickaH"} See help for \code{MVLcc()}
#' }

#' @examples 
#' xi <- heather$coarse
#' xiim <- as.im(xi, value = TRUE, na.replace = FALSE)
#' mvlests <- mvl(xiim, seq(1, 10, by = 0.1))

mvl <- function(xiim, boxwidths,
                           estimators = c("MVLg.mattfeldt", "MVLg.pickaint", "MVLg.pickaH",
                                          "MVLcc.mattfeldt", "MVLcc.pickaint",
                                          "MVLc", "MVLgb"),
                includenormed = FALSE,
                includepaircorr = FALSE,
                includecovar = FALSE){
  mvlgestimaterequests <- estimators %in% MVLgestimatornames
  mvlccestimaterequests <- estimators %in% MVLccestimatornames
  
  phat <- coverageprob(xiim)
  cvchat <- racscovariance(xiim, setcov_boundarythresh = 0.1 * area.owin(solutionset(is.na(xiim))))
  mvlgs <- mvlccs <- mvlc.est <- mvlgb.est <- NULL
  cpp1 <- NULL
  if (sum(mvlgestimaterequests) + sum(mvlccestimaterequests) > 0){
    cpp1 <- cppicka(xiim, setcov_boundarythresh = 0.1 * area.owin(solutionset(is.na(xiim))))
  }
  #function that computes the covariance-based estimates of MVL
  mvlcovarbased <- mvl.cvchat(boxwidths = boxwidths, estimators = estimators, phat = phat, cvchat = cvchat, cpp1 = cpp1)
  
  #the MVLgb estimate
  if ("MVLgb" %in% estimators){
    mvlgb.est <- mvlgb(sidelengths = boxwidths, xiim = xiim)
    if (sum(!vapply(mvlgb.est[,fvnames(mvlgb.est), drop = TRUE], is.na, FUN.VALUE = TRUE)) < 2){
      warning("mvlgb() returns estimates for 1 or fewer of the provided box widths. Results from mvlgb() will be ignored from the final results.")
      mvlgb.est <- NULL
    }
  }
  
  mvl.ests <- c(mvlcovarbased, list(mvlgb = mvlgb.est))
  mvl.ests <- mvl.ests[!vapply(mvl.ests, is.null, FUN.VALUE = FALSE)]
  if (any(!vapply(mvl.ests[-1], function(x) compatible.fv(A = mvl.ests[[1]], B = x), FUN.VALUE = FALSE))){
    warning("Some MVL estimates have differing argument values. These will be harmonised.")
    mvl.ests <- harmonise.fv(mvl.ests)
  }
  mvls.fv <- collapse.fv(mvl.ests, different = "MVL")
  names(mvls.fv) <- c(fvnames(mvls.fv, ".x"), names(mvl.ests))
  
  allfvs <- list(mvlests = mvls.fv)
  
  if (includenormed){
    #compute MVLs normalised at zero
    normdmvls <- eval.fv(mvls.fv / ( phat * (1 - phat) / phat^2), relabel = FALSE)
    normdmvls <- prefixfv(normdmvls,
                     tagprefix="n_",
                     descprefix="normalised at zero",
                     lablprefix="plain(nrmd)~")
    allfvs <- c(allfvs, list(normdmvls = normdmvls))
  }
  
  if (includecovar || includepaircorr){
    #computing an isotropic covariance
    if (sum(mvlgestimaterequests) + sum(mvlccestimaterequests) > 0){
      impcovar <- balancedracscovariance.cvchat(cvchat, cpp1 = cpp1, phat = phat, modification = "pickaH")
      isocovar <- rotmean(impcovar, padzero = FALSE, Xname = "covar", result = "fv")
    } else {
      isocovar <- rotmean(cvchat, padzero = FALSE, Xname = "covar", result = "fv")
    }
    isocovar <- tweak.fv.entry(isocovar, "f", new.labl = "C(r)", new.desc = "isotropic covariance", new.tag = "C")
    isocovar <- rebadge.fv(isocovar,
                           new.ylab = "Isotropic Covariance",
                           new.fname = "C(r)")
    allfvs <- c(allfvs, list(covar = isocovar))
  }
  
  if (includepaircorr){
    #compute isotropic pair correlation
    if (sum(mvlgestimaterequests) + sum(mvlccestimaterequests) > 0){
      pclnest <- pclns.cvchat(cvchat, cpp1 = cpp1, phat = phat, modifications = "pickaH")[[1]]
      isopcln <- rotmean(pclnest, padzero = FALSE, Xname = "pcln", result = "fv")
    } else {
      isopcln <- eval.fv(isocovar / phat^2, relabel = TRUE) #if no improvements available then use traditional estimates
  
    }
    isopcln <- tweak.fv.entry(isopcln, "f", new.labl = "g(r)", new.desc = "isotropic pair-correlation", new.tag = "g")
    isopcln <- rebadge.fv(isopcln,
                          new.ylab = "Pair-Correlation",
                          new.fname = "g(r)")
    allfvs <- c(allfvs, list(paircorr = isopcln))
  }

  return(allfvs)
}

mvl.cvchat <- function(boxwidths,
                      estimators = c("MVLg.mattfeldt", "MVLg.pickaint", "MVLg.pickaH",
                                     "MVLcc.mattfeldt", "MVLcc.pickaint",
                                     "MVLc", "MVLgb"),
                      phat = NULL,
                      cvchat = NULL,
                      cpp1 = NULL){
  mvlgestimaterequests <- estimators %in% MVLgestimatornames
  mvlccestimaterequests <- estimators %in% MVLccestimatornames
  mvlgs <- mvlccs <- mvlc.est <- mvlgb.est <- NULL
  
  if (sum(mvlgestimaterequests) + sum(mvlccestimaterequests) > 0){
    stopifnot(!is.null(cpp1), !is.null(cvchat), !is.null(phat))
    if (sum(mvlgestimaterequests) > 0){
      pcln.ests <- pclns.cvchat(cvchat, cpp1 = cpp1, phat = phat, modifications = gsub("MVLg.", "", estimators[mvlgestimaterequests]))
      mvlgs <- lapply(pcln.ests, FUN = mvlg, boxes = boxwidths)
    }
    if (sum(mvlccestimaterequests) > 0){
      ccvc.ests <- ccvcs.cvchat(cvchat, cpp1, phat, modifications = gsub("MVLcc.", "", estimators[mvlccestimaterequests]))
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
