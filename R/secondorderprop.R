#' @title Estimate Second-Order Properties of a RACS
#' @description Estimates many second order properties of RACS, mass variance lacunarity, covariance, centred covariance,#' and pair-correlation.
#' This can be faster than computing estimates of multiple second order properties separately as  
#' Fourier transforms of the binary map are not repeated.

#' @export secondorderprops
#' @param xiim A \pkg{spatstat} \code{im} object with pixel values that are either TRUE, FALSE or NA. TRUE represents foreground, FALSE respresents background and NA represents unobserved locations.
#' @param mvlargs Arguments passed to \code{mvlgb} and \code{mvl.cvchat}. If NULL then MVL will not be estimated.
#' You can also request to 
#' @param covarargs Arguments passed to \code{racscovariance.cvchat}. If NULL then covariance will not be returned.
#' @param cencovarargs NOT YET IMPLEMENTED
#' @param paircorrargs Arguments passed to \code{paircorr.cvchat}. If NULL then pair correlation will not be returned.
#' @param returnrotmean Logical. If FALSE the anisotropic estimates of covariance and pair-correlation will be returned as \code{im} objects.
#' If TRUE then average covariance and pair-correlation over all directions will be returned as \code{fv} objects.

#' @examples 
#' xi <- heather$coarse
#' xiim <- as.im(xi, value = TRUE, na.replace = FALSE)
#' mvlargs = list(boxwidths = seq(1, 10, by = 0.1), estimators = c("MVLgb", "MVLc"))
#' covarargs = list(modifications = "all")
#' paircorrargs = list(modifications = "pickaH")
#' returnrotmean = TRUE
#' secondests <- secondorderprops(xiim, mvlargs = mvlargs, covarargs = covarargs, paircorrargs = paircorrargs, returnrotmean = FALSE)
#'#to test, with and without returnrotmean, with various parts NULL, with only MVLgb, only MVLc.pickaH or no estimators.
#'#with discs

secondorderprops <- function(xiim,
                             mvlargs = NULL,
                             covarargs = NULL,
                             cencovarargs = NULL,
                             paircorrargs = NULL,
                             returnrotmean = FALSE
                            ){
  cvchatT <- tradcovarest(xiim)
  cpp1 <- cppicka(xiim)
  phat <- coverageprob(xiim)
  
  outlist <- list()
  
  #first compute MVL ests copying MVL()
  if (!is.null(mvlargs)){ 
    mvl.ests <- list()
    if  (any("MVLgb" != mvlargs[["estimators"]])){
      mvlcovarbased <- do.call(mvl.cvchat, args = c(mvlargs, list(phat = phat, cvchat = cvchatT, cpp1 = cpp1)))
      mvl.ests <- c(mvl.ests, mvlcovarbased)
    }
    if ("MVLgb" %in% mvlargs[["estimators"]]) {
      mvlgb.est <- mvlgb(sidelengths = mvlargs[["boxwidths"]], xiim = xiim)
      if (sum(!vapply(mvlgb.est[,fvnames(mvlgb.est), drop = TRUE], is.na, FUN.VALUE = TRUE)) < 2){
        warning("mvlgb() returns estimates for 1 or fewer of the provided box widths. Results from mvlgb() will be ignored from the final results.")
        mvlgb.est <- NULL
      }
      mvl.ests <- c(mvl.ests, list(mvlgb = mvlgb.est))
    }
    #combind the mvl ests
    mvl.ests <- mvl.ests[!vapply(mvl.ests, is.null, FUN.VALUE = FALSE)]
    if (any(!vapply(mvl.ests[-1], function(x) compatible.fv(A = mvl.ests[[1]], B = x), FUN.VALUE = FALSE))){
      warning("Some MVL estimates have differing argument values. These will be harmonised.")
      mvl.ests <- harmonise.fv(mvl.ests)
    }
    mvls.fv <- collapse.fv(mvl.ests, different = "MVL")
    names(mvls.fv) <- c(fvnames(mvls.fv, ".x"), names(mvl.ests))
    outlist <- c(outlist, list(MVL = mvls.fv))
  }
  ### mvl estimation finished ###
  
  #covariance computations
  if (!is.null(covarargs)) {
    cvchats <- do.call(racscovariance.cvchat, args = c(list(cvchat = cvchatT, cpp1 = cpp1, phat = phat), covarargs, drop = FALSE))
    if (returnrotmean){
      isocovars <- lapply(cvchats, rotmean, padzero = FALSE, Xname = "covar", result = "fv")
      isocovars <- lapply(isocovars, function(x) {
        x <- tweak.fv.entry(x, "f", new.labl = "C(r)", new.desc = "isotropic covariance", new.tag = "C")
        return(x)
      })
      isocovars <- collapse.fv(isocovars, different = "C")
      cvchats <- isocovars
    }
    outlist <- c(outlist, list(covariance = cvchats))
  }
  
  if (!is.null(paircorrargs)){
    #compute isotropic pair correlation
    pclnests <- do.call(paircorr.cvchat, c(list(cvchat = cvchatT, cpp1 = cpp1, phat = phat), paircorrargs, drop = FALSE))
    if (returnrotmean){
      isopclns <- lapply(plcnests, rotmean, padzero = FALSE, Xname = "covar", result = "fv")
      isopclns <- lapply(isopclns, function(x) {
        x <- tweak.fv.entry(x, "f", new.labl = "g(r)", new.desc = "isotropic pair-correlation", new.tag = "g")
        return(x)
      })
      isopclns <- collapse.fv(isopclns, different = "g")
      pclnests <- isopclns
    }
    outlist <- c(outlist, list(paircorr = pclnests))
  }
  
  return(outlist)
}
                             
