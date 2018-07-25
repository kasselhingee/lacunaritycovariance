#' @title Median Integrated Squared Error of a List of fv Objects
#' @export median_ise.fvlist
#' @description
#' This function assumes that fv objects are each realisations of the same (stochastic) object. 
#' It returns pointwise summaries such as observed sample mean and sample variance.

#' @author Kassel Hingee

#' @param object List of function objects to compare to benchfv
#' @param benchfv An fv function object to treat as the true or reference values.
#' @param domainlim A list of two entries giving the lower and upper bound of the domain to integrate over
#' @param equiv Passed to eval.fv - see help for eval.fv for more details.
#'        It is gives the functional value names to map to each other when comparing to benchfv. Default is NULL.
#' @param acceptableISEerrorrate The number of bad integrations that it is ok to ignore when computing the median ISE.
#' @param ... Arguments passed to integrate().
#' @return An fv object containing the pointwise mean, pointwise sample variance, pointwise maximum and pointwise minimum.

#' @examples
#' obspatterns <- rpoispp(10, nsim = 10)
#' object <- lapply(obspatterns, Kest, W = Frame(obspatterns[[1]]))
#' benchfv <- as.fv(object[[1]][, c("r", "theo"), drop = TRUE])
#' names(benchfv) <- c("r", "bench")
#' median_ise.fvlist(object, benchfv, c(0.1, 0.2), equiv = list(bench = "iso"))
#' 
#' @details 
#' The function fails if there is y-value name of the reference fv object
#'  is equal to a y-value name in the list fv objects that that you don't want to compare to
#'   (e.g. if the list of fv objects also contain the reference value).
median_ise.fvlist <- function(object, benchfv, domainlim, equiv = NULL, acceptableISEerrorrate = 0.1, ...){
  isel <- lapply(object, ise, benchfv = benchfv, domainlim = domainlim, equiv = equiv, ...)
  ermessageok <- grepl("^OK$", lapply(isel, "[[", "message"))
  if (is.null(acceptableISEerrorrate)){ acceptableISEerrorrate <- 0.1 }
  if ( sum(ermessageok) / length(isel) < 1 - acceptableISEerrorrate){
    #warning(lapply(isel[!ermessageok], "[[", "message"))
    warning(sprintf("The fraction of integrated squared error calls that had a warning was %.2g. This is greater than acceptableISEerrorrate = %3.2g.", 1 - sum(ermessageok) / length(isel),  acceptableISEerrorrate))
    return(NA)
  }
  iselvalues <- vapply(isel, "[[", "value", FUN.VALUE = 0.0)
  return( median(iselvalues, na.rm = TRUE) )
}



ise <- function(fvobj, benchfv, domainlim, equiv = NULL, ...){
  harmfvs <- harmonise.fv(fvobj, benchfv)
  sefv <- eval.fv( (a - b)^2, envir = list(a = harmfvs[[1]], b = harmfvs[[2]]), equiv = equiv, relabel = FALSE)
  se.fun <- as.function.fv(sefv)
  #following code makes the integration
  # between known values for both functions
  finitevals <-
    is.finite(sefv[, fvnames(sefv, ".y"), drop = TRUE])
  if (min(sefv[, fvnames(sefv, ".x"), drop = TRUE][finitevals]) > domainlim[[1]]){
    return(list(value = NA, message = "integration domain larger the function domain"))
  }
  if (max(sefv[, fvnames(sefv, ".x"), drop = TRUE][finitevals]) < domainlim[[2]]){
    return(list(value = NA, message = "integration domain larger the function domain"))
  }
  lower <- max( min(sefv[, fvnames(sefv, ".x"), drop = TRUE][finitevals]), domainlim[[1]])
  upper <- min( max(sefv[, fvnames(sefv, ".x"), drop = TRUE][finitevals]), domainlim[[2]])
  intse <- integrate(se.fun, lower, upper, stop.on.error = FALSE, ...)
  #subdivisions error from MVL going too high too fast close to 0?! I have ignored these curves by setting them to NA if there are few enough of them.
  # stopifnot( (intse$message == "maximum number of subdivisions reached")
  #            || (intse$message == "extremely bad integrand behaviour")
  #            || (intse$message == "OK"))
  if (intse$message != "OK"){ intse$value <- NA}
  intse$value <- intse$value #/ (upper - lower)
  intse$abs.error <-  intse$abs.error #/ (upper - lower)
  return(intse)
}
