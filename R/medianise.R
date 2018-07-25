#' @title Median and Mean Integrated Squared Error of a List of fv Objects
#' @export median_ise.fvlist  mean_ise.fvlist
#' @description
#' Computes the median integrated squared error (defined below) of a list of fv objects given reference fv object.

#' @author Kassel Hingee

#' @param object List of function objects to compare to reffv
#' @param reffv An fv function object to treat as the true or reference values.
#' @param domainlim A list of two entries giving the lower and upper bound of the domain to integrate over
#' @param equiv Passed to eval.fv - see help for eval.fv for more details.
#'        It is gives the functional value names to map to each other when comparing to reffv. Default is NULL.
#' @param acceptableISEerrorrate The number of bad integrations that it is ok to ignore when computing the median ISE.
#' @param avoverdomain  Divide the integral by the length of the domain (i.e. divide by \eqn{b_i} - \eqn{a_i} for each \eqn{\hat{f}_i}
#' @param fixeddomain If the domain of the function \eqn{\hat{f}_i}
#'  is a non-equal subset of the domain given by \code{domainlim} then set the integrated error to NA.
#' @param ... Arguments passed to integrate().
#' @return A numeric value that is the Median Integrated Squared Error - see Details.

#' @examples
#' obspatterns <- rpoispp(10, nsim = 10)
#' object <- lapply(obspatterns, Kest, W = Frame(obspatterns[[1]]))
#' reffv <- as.fv(object[[1]][, c("r", "theo"), drop = TRUE])
#' names(reffv) <- c("r", "ref")
#' median_ise.fvlist(object, reffv, c(0.1, 0.2), equiv = list(ref = "iso"))
#' mean_ise.fvlist(object, reffv, c(0.1, 0.2), equiv = list(ref = "iso"))
#' 
#' @details 
#' #' We define the integrated squared error of a function \eqn{\hat{f}}
#' and a reference function \eqn{f} as 
#' \deqn{
#' ISE(\hat{f}) := \int_{a_i}^{b_i} \left(\hat{f}(s) - f(s) \right)^2 ds,
#' }
#' where the domain of integration $[a_i, b_i]$ is 
#' the largest possible interval between \code{domainlim[[1]]} and \code{domainlim[[2]]} for which the function
#' \eqn{\hat{f}} exists.
#' The integral is computed using the numerical integration function \code{integrate} and
#'  if this numerical integration fails then ISE(\hat{f}) is given a value of NA.
#' (Unless \code{fixeddomain} is \code{TRUE} in which case the domain of integration is fixed to \code{domainlim})
#' 
#' \code{median_ise.fvlist} and \code{meain_ise.fvlist} computes the median and mean, respectively, of this integrated squared error for a collection
#' of functions \eqn{\hat{f}_i}, \eqn{i = 1, 2, 3, ... n}.
#' 
#' \code{acceptableISEerrorrate} can be used to set the proportion of ISE that are NA valued
#'  that is acceptable for calculating the median or mean.
#'  The mean and median will both ignore NA values.
#' 
#' The function fails if there is y-value name of the reference fv object
#'  is equal to a y-value name in the list fv objects that that you don't want to compare to
#'   (e.g. if the an fv object in the list also contains the reference value).
median_ise.fvlist <- function(object, reffv, domainlim, equiv = NULL,
                              avoverdomain = FALSE, fixeddomain = FALSE, acceptableISEerrorrate = 0.1, ...){
  isel <- lapply(object, ise, reffv = reffv, domainlim = domainlim, equiv = equiv,
                 avoverdomain = avoverdomain, fixeddomain = fixeddomain, ...)
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

#' @describeIn median_ise.fvlist The mean integrated squared error.
mean_ise.fvlist <- function(object, reffv, domainlim, equiv = NULL,
                              avoverdomain = FALSE, fixeddomain = FALSE, acceptableISEerrorrate = 0, ...){
  isel <- lapply(object, ise, reffv = reffv, domainlim = domainlim, equiv = equiv,
                 avoverdomain = avoverdomain, fixeddomain = fixeddomain, ...)
  ermessageok <- grepl("^OK$", lapply(isel, "[[", "message"))
  if (is.null(acceptableISEerrorrate)){ acceptableISEerrorrate <- 0 }
  if ( sum(ermessageok) / length(isel) < 1 - acceptableISEerrorrate){
    #warning(lapply(isel[!ermessageok], "[[", "message"))
    warning(sprintf("The fraction of integrated squared error calls that had a warning was %.2g. This is greater than acceptableISEerrorrate = %3.2g.", 1 - sum(ermessageok) / length(isel),  acceptableISEerrorrate))
    return(NA)
  }
  iselvalues <- vapply(isel, "[[", "value", FUN.VALUE = 0.0)
  return( mean(iselvalues, na.rm = TRUE) )
}

ise <- function(fvobj, reffv, domainlim, equiv = NULL, avoverdomain = FALSE, fixeddomain = FALSE, ...){
  harmfvs <- harmonise.fv(fvobj, reffv)
  withCallingHandlers(
    sefv <- eval.fv( (a - b)^2, envir = list(a = harmfvs[[1]], b = harmfvs[[2]]), equiv = equiv, relabel = FALSE),
      warning = function(w){
        if(grepl("enforcing compatibility", w$message)){
          invokeRestart( "muffleWarning" )
        } else {
        }
      })
  
  se.fun <- as.function.fv(sefv)
  finitevals <- is.finite(sefv[, fvnames(sefv, ".y"), drop = TRUE])
  if (fixeddomain){#if domain is fixed then throw out NAs when the function objects aren't long enough
    if (min(sefv[, fvnames(sefv, ".x"), drop = TRUE][finitevals]) > domainlim[[1]]){
      return(list(value = NA, message = "integration domain larger than function domain"))
    }
  }
  lower <- max( min(sefv[, fvnames(sefv, ".x"), drop = TRUE][finitevals]), domainlim[[1]])
  upper <- min( max(sefv[, fvnames(sefv, ".x"), drop = TRUE][finitevals]), domainlim[[2]])
  intse <- integrate(se.fun, lower, upper, stop.on.error = FALSE, ...)
  if (intse$message != "OK"){ intse$value <- NA}
  if (avoverdomain) { intse$value <- intse$value / (upper - lower) }
  if (avoverdomain) { intse$abs.error <-  intse$abs.error / (upper - lower) }
  return(intse)
}
