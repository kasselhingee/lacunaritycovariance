#' @title Estimate Second-Order Properties of a RACS
#' @description Estimates many second order properties of RACS, gliding box lacunarity, covariance, centred covariance,
#' and pair-correlation.
#' This can be faster than computing estimates of multiple second order properties separately as  
#' Fourier transforms of the binary map are not repeated.

#' @export secondorderprops
#' @param xiim A \pkg{spatstat} \code{im} object with pixel values that are either TRUE, FALSE or NA. TRUE represents foreground, FALSE represents background and NA represents unobserved locations.
#' @param gblargs Arguments passed to \link{\code{gblemp}} and \link{\code{gbl.cvchat}}. If NULL then GBL will not be estimated.
#' You can also request to 
#' @param covarargs Arguments passed to \code{racscovariance.cvchat}. If NULL then covariance will not be returned.
#' @param cencovarargs NOT YET IMPLEMENTED
#' @param paircorrargs Arguments passed to \code{paircorr.cvchat}. If NULL then pair correlation will not be returned.
#' @param returnrotmean Logical. If FALSE the anisotropic estimates of covariance and pair-correlation will be returned as \code{im} objects.
#' If TRUE then average covariance and pair-correlation over all directions will be returned as \code{fv} objects.

#' @examples 
#' xi <- heather$coarse
#' xiim <- as.im(xi, value = TRUE, na.replace = FALSE)
#' gblargs = list(boxwidths = seq(1, 10, by = 1), estimators = c("GBLgb", "GBLc"))
#' covarargs = list(estimators = "all")
#' paircorrargs = list(estimators = "pickaH")
#' returnrotmean = TRUE
#' secondests <- secondorderprops(xiim,
#'    gblargs = gblargs,
#'    covarargs = covarargs,
#'    paircorrargs = paircorrargs, 
#'    returnrotmean = FALSE)

secondorderprops <- function(xiim,
                             gblargs = NULL,
                             covarargs = NULL,
                             cencovarargs = NULL,
                             paircorrargs = NULL,
                             returnrotmean = FALSE
                            ){
  cvchatT <- plugincvc(xiim)
  cpp1 <- cppicka(xiim)
  phat <- coverageprob(xiim)
  
  outlist <- list()
  
  #first compute GBL ests copying GBL()
  if (!is.null(gblargs)){
    gbl.ests <- list()
    if  (any("GBLgb" != gblargs[["estimators"]]) || is.null(gblargs[["estimators"]]) || ("all" %in% gblargs[["estimators"]])){
      gblcovarbased <- do.call(gbl.cvchat, args = c(gblargs, list(phat = phat, cvchat = cvchatT, cpp1 = cpp1)))
      gbl.ests <- c(gbl.ests, gblcovarbased)
    }
    if (("GBLgb" %in% gblargs[["estimators"]]) || is.null(gblargs[["estimators"]]) || ("all" %in% gblargs[["estimators"]])) {
      gblemp.est <- gblemp(boxwidths = gblargs[["boxwidths"]], xiim = xiim)
      if (sum(!vapply(gblemp.est[,fvnames(gblemp.est), drop = TRUE], is.na, FUN.VALUE = TRUE)) < 2){
        warning("gblemp() returns estimates for 1 or fewer of the provided box widths. Results from gblemp() will be ignored from the final results.")
        gblemp.est <- NULL
      }
      gbl.ests <- c(gbl.ests, list(gblemp = gblemp.est))
    }
    #combind the gbl ests
    gbl.ests <- gbl.ests[!vapply(gbl.ests, is.null, FUN.VALUE = FALSE)]
    if (is.owin(gblargs[["boxwidths"]][[1]])){
      gbls <- do.call(cbind, args = gbl.ests)
    } else {
      if (any(!vapply(gbl.ests[-1], function(x) compatible.fv(A = gbl.ests[[1]], B = x), FUN.VALUE = FALSE))){
        warning("Some GBL estimates have differing argument values. These will be harmonised.")
        gbl.ests <- harmonise.fv(gbl.ests)
      }
      gbls <- collapse.fv(gbl.ests, different = "GBL")
      names(gbls) <- c(fvnames(gbls, ".x"), names(gbl.ests))
    }
    outlist <- c(outlist, list(GBL = gbls))
  }
  ### gbl estimation finished ###
  
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
      isopclns <- lapply(pclnests, rotmean, padzero = FALSE, Xname = "covar", result = "fv")
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
                             
