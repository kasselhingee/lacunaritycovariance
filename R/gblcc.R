#' @title Centred covariance based estimates of gliding box lacunarity
#' @export gblcc gblcc.inputcovar
#'
#' @description Estimates the gliding box lacunarity (GBL) of a stationary RACS using centred covariance estimates (Hingee et al., 2017).
#'  The centred covariance and coverage probability can be provided or estimated from binary map.

#' @details If we denote the estimated centred covariance by
#' \eqn{\hat{\kappa}(v)}{k(v)} and coverage probability \eqn{\hat{p}}{p} then
#' the estimate of GBL is
#' \deqn{1 + \frac{1}{\hat{p}^2 |B|^2}\int \gamma_B(v)\hat{\kappa}(v)dv, }{1 + \int gammaB(v) k(v) dv  /  (p^2 |B|^2),}
#' where \eqn{B} is each of the sets (often called a box) specified by \code{boxes},
#' \eqn{\gamma_B}{gammaB} is the set covariance of \eqn{B},
#' \eqn{|B|} is the area of \eqn{B},
#' \eqn{p} is the coverage probability of a stationary RACS.
#' 
#' The set covariance of \eqn{B} is computed empirically using \pkg{spatstat}'s \code{\link[spatstat.geom]{setcov}} function, which converts \eqn{B} into a binary pixel mask using \code{\link[spatstat.geom]{as.mask}} defaults. Computation speed can be increased by setting a small default number of pixels, \code{npixel}, in \pkg{spatstat}'s global options (accessed through \code{\link[spatstat.geom]{spatstat.options}}), however fewer pixels also decreases the accuracy of the GBL computation.


#' @param boxes Either a list of side lengths for square boxes or a list of \code{owin} objects of any shape.
#' @param cencovar  A \code{im} object containing the centred covariance function
#' @param p The coverage probability. Typically estimated by the fraction of the observation window covered by the set of interest.
#' @param xiim An observation of a stationary RACS as an \code{im} object. \code{xiim} must have values of either 1, 0 or NA; 1 denotes inside the RACS, 0 denotes outside, and NA denotes unobserved.
#' @param estimator If an observation \code{xiim} is passed then \code{estimator} will select the balancing method that \code{ccvc} uses to estimate the centred covariance.
#' @param integrationMethod The integration method used by [innerprod.im()]. Default is 'harmonisesum' due centred covariance approaching zero for large vectors.

#' @return If \code{boxes} is a list of numerical values then GBL is estimated for square boxes with side length given by \code{boxes}.
#'  The returned object is then an \code{fv} object containing estimates of GBL, box mass variance and box mass mean.
#'  
#'  If \code{boxes} is a list of \code{owin} objects then \code{gblcc} returns a dataframe of with columns corresponding to estimates of GBL, box mass variance and box mass mean.
#'  Note if \code{NA} or \code{NaN} values in the \code{covariance} object are used then \code{gblc} will return \code{NA} or \code{NaN} instead of an GBL value. 

#' @references
#' Hingee K, Baddeley A, Caccetta P, Nair G (2019). Computation of lacunarity from covariance of spatial binary maps. \emph{Journal of Agricultural, Biological and Environmental Statistics}, 24, 264-288. DOI: 10.1007/s13253-019-00351-9.

#' @examples
#' xi <- heather$coarse
#' cencovar <- cencovariance(xi, obswin = Frame(xi), estimators = c("pickaH"), drop = TRUE)
#' p <- area(xi) / area(Frame(xi))
#' sidelengths <- seq(0.3, 14, by = 1)
#' 
#' # reduce resolution in setcov() for faster (less accurate) computation 
#' oldopt <- spatstat.options()
#' spatstat.options("npixel" = 2^5)
#'
#' # compute GBL estimate for square boxes from estimated centred covariance
#' gblccest <- gblcc(sidelengths, cencovar, p)
#' 
#' # compute GBL estimate for boxes that are discs
#' discboxes <- lapply(sidelengths / 2, disc)
#' discgbls <- gblcc(discboxes, cencovar, p)
#' 
#' # compute GBL estimates from binary map
#' xiim <- as.im(xi, na.replace = 0)
#' gblccest <- gblcc(sidelengths, xiim = xiim, estimator = "pickaH")
#'
#' spatstat.options(oldopt)
#' 
#' @keywords spatial nonparametric 
gblcc <- function(boxes, cencovar = NULL, p = NULL, xiim = NULL, estimator = "pickaH", integrationMethod = "harmonisesum"){
  if (!(is.null(cencovar) && is.null(p))){
    if (!is.null(xiim)){stop("xiim (an observation image) and cencovar or p were given. Either cencovar and p must be supplied or xiim supplied.")}
    lacv <- gblcc.inputcovar(boxes, cencovar, p)
    unitname <- unitname(cencovar)
  } else if (!is.null(xiim)){
    p <- sum(xiim) / sum(is.finite(xiim$v))
    cencovar <- cencovariance(xiim, estimators = c(estimator), drop = TRUE)
    lacv <- gblcc.inputcovar(boxes, cencovar, p, integrationMethod = integrationMethod)
    unitname <- unitname(xiim)
  } else {
    stop("Input requires specification of xiim or cencovar and p")
  }

  lacsdf <- as.data.frame(lacv)
  if (mode(boxes) %in% c("integer", "numeric")){
    lacsdf <- cbind(s = boxes, lacsdf)
    #recommended xlim:
    alim.min <- 1
    alim.max <- min(which(vapply(lacsdf[, "GBL"], is.na, FUN.VALUE = TRUE)), nrow(lacsdf))
    lacfv <- fv(lacsdf,
                argu = "s",
                valu = "GBL",
                fmla = ".y ~ s",
                alim = c(lacsdf[alim.min, "s"], lacsdf[alim.max, "s"]),
                ylab = quote(GBL[c]),
                unitname = unitname(xiim),
                labl = c("Box Width",
                         "GBL",
                         "Var(BoxMass)",
                         "Mean[BoxMass]"),
                desc = c("side lengths of boxes", 
                         "GBL derived from covariance",
                         "A covariance-based estimate of the variance in box mass",
                         "An estimate of mean box mass using coverage probability"
                ),
                fname = "GBL"
    )
    fvnames(lacfv, a = ".") <- "GBL"
    return(lacfv)
  } else (return(lacsdf))
}

#' @describeIn gblcc GBL estimates from already estimated centred covariance.
gblcc.inputcovar <- function(boxes, cencovar, p, integrationMethod =  "harmonisesum"){
  stopifnot(is.im(cencovar))
  stopifnot(is.numeric(p))
  if (mode(boxes) %in% c("integer", "numeric")){
    squares <- lapply(boxes, square) #make into owin rectangles
    boxcov <- lapply(squares, setcov) #setcov is analytic for squares according to help, couldn't see it in code though.
                                     #regardless - it produces much better plots then my own function setcovsquare did
    boxarea <- boxes ^ 2
  }
  else { #box must be a list of owin objects
    boxcov <- lapply(boxes, setcov) #numerical
    boxarea <- lapply(boxes, area.owin)
    boxarea <- unlist(boxarea)
  }

  integrationresults <- mapply(innerprod.im, boxcov, list(cencovar), outsideA = 0, outsideB = NA, na.replace = TRUE, method = integrationMethod, SIMPLIFY = FALSE) # the list around the cencovar is necessary to stop mapply unlisting the image itself

  coefvar2 <- unlist(integrationresults) / (p ^ 2 * boxarea ^ 2)
  return(list(
    GBL = coefvar2 + 1,
    s2 = unlist(integrationresults),
    xbar = p * boxarea
  ))
}
