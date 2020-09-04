#' @title Estimate Second-Order Properties of a RACS
#' @description Estimates many second order properties of RACS, gliding box lacunarity, covariance, centred covariance,
#' and pair-correlation.
#' This can be faster than computing estimates of multiple second order properties separately as  
#' Fourier transforms of the binary map are not repeated.

#' @export secondorderprops
#' @param xiim A \pkg{spatstat} \code{im} object with pixel values that are either TRUE, FALSE or NA. TRUE represents foreground, FALSE represents background and NA represents unobserved locations.
#' @param gblargs A list of named arguments passed to \code{\link{gblemp}} and \code{\link{gbl.cvchat}}. The estimators used can be specified by passing an argument named 'estimators' contain a list of estimator names, as given in \code{\link{gbl}}. If NULL then GBL will not be estimated.
#' @param covarargs A list of named arguments passed to \code{\link{racscovariance.cvchat}}. If NULL then covariance will not be returned.
#' @param cencovarargs A list of named arguments passed to \code{\link{cencovariance.cvchat}}. If NULL then pair correlation will not be returned.
#' @param paircorrargs A list of named arguments passed to \code{\link{paircorr.cvchat}}. If NULL then pair correlation will not be returned.
#' @param returnrotmean Logical. If FALSE the anisotropic estimates of covariance and pair-correlation will be returned as \code{im} objects.
#' If TRUE then average covariance and pair-correlation over all directions will be returned as \code{fv} objects.

#' @examples 
#' xiim <- as.im(heather$coarse, value = TRUE,
#'               na.replace = FALSE)
#' gblargs = list(boxwidths = seq(1, 10, by = 1), estimators = c("GBLemp", "GBLcc.pickaH"))
#' covarargs = list(estimators = "all")
#' cencovarargs = list(estimators = "pickaH")
#' paircorrargs = list(estimators = "pickaH")
#' returnrotmean = TRUE
#' secondests <- secondorderprops(xiim,
#'    gblargs = gblargs,
#'    covarargs = covarargs,
#'    cencovarargs = cencovarargs,
#'    paircorrargs = paircorrargs, 
#'    returnrotmean = FALSE)

#' @section Warning: The user interface for this function is substantially more  stretched the knowledge of the author, Kassel Hingee. Therefore, there is greater chance of encountering bugs. Kassel Hingee apologises for any bugs you encounter, and requests to be informed (thank you!).

#' @return 
#' A named list of estimated properties. When multiple estimators have been requested for the same property, then the entry in the list is itself a list, with each entry corresponding to a different estimator.


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
    if  (any("GBLemp" != gblargs[["estimators"]]) || is.null(gblargs[["estimators"]]) || ("all" %in% gblargs[["estimators"]])){
      gblcovarbased <- do.call(gbl.cvchat, args = c(gblargs, list(phat = phat, cvchat = cvchatT, cpp1 = cpp1)))
      gbl.ests <- c(gbl.ests, gblcovarbased)
    }
    if (("GBLemp" %in% gblargs[["estimators"]]) || is.null(gblargs[["estimators"]]) || ("all" %in% gblargs[["estimators"]])) {
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
      isocovars <- collapse.fv(harmonise.fv(isocovars), different = "C")
      cvchats <- isocovars
    }
    outlist <- c(outlist, list(covariance = cvchats))
  }
  
#centred covariance computations
  if (!is.null(cencovarargs)) {
    ccvchats <- do.call(cencovariance.cvchat, args = c(list(cvchat = cvchatT, cpp1 = cpp1, phat = phat), cencovarargs, drop = FALSE))
    if (returnrotmean){
      isocencovars <- lapply(ccvchats, rotmean, padzero = FALSE, Xname = "cencovar", result = "fv")
      isocencovars <- lapply(isocencovars, function(x) {
        x <- tweak.fv.entry(x, "f", new.labl = "k(r)", new.desc = "isotropic centred covariance", new.tag = "C")
        return(x)
      })
      isocencovars <- collapse.fv(harmonise.fv(isocencovars), different = "C")
      ccvchats <- isocencovars
    }
    outlist <- c(outlist, list(cencovariance = ccvchats))
  }
  
  if (!is.null(paircorrargs)){
    pclnests <- do.call(paircorr.cvchat, c(list(cvchat = cvchatT, cpp1 = cpp1, phat = phat), paircorrargs, drop = FALSE))
    #compute isotropic pair correlation
    if (returnrotmean){
      isopclns <- lapply(pclnests, rotmean, padzero = FALSE, Xname = "paircorr", result = "fv")
      isopclns <- lapply(isopclns, function(x) {
        x <- tweak.fv.entry(x, "f", new.labl = "g(r)", new.desc = "isotropic pair-correlation", new.tag = "g")
        return(x)
      })
      isopclns <- collapse.fv(harmonise.fv(isopclns), different = "g")
      pclnests <- isopclns
    }
    outlist <- c(outlist, list(paircorr = pclnests))
  }
  
  return(outlist)
}
                             
