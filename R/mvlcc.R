#' @title Centred covariance based calculations of mass variance lacunarity
#' @export mvlcc
#'
#' @description Estimates the mass variance lacunarity (MVL) of a stationary RACS from a bi-tonal image using centred covariance estimates.
#'  The centred covariance and coverage probability can be provided or estimated from an image.

#' @details
#' \code{mvl} and \code{lac} are aliases of mvlc.
#' If we denoted the estimated centred covariance by \eqn{\hat{\kappa}(v)} and coverage probability \eqn{\hat{p}} then the estimate of MVL is
#' \deqn{\frac{1}{\hat{p}^2 |B|^2}\int \gamma_B(v)\hat{kapp}(v)dv}

#' @param boxes Either a list of sidelengths for square boxes or a list of \code{owin} objects of any shape.
#' @param cencovar  A \code{im} object containing the centred covariance function
#' @param p The coverage probability. Typically estimated by the fraction of the observation window covered by the set of interest.
#' @param xiim An observation of a stationary RACS in \code{im} format. \code{xiim} must have values of either 1, 0 or NA; 1 denotes inside the RACS, 0 denotes outside, and NA denotes unobserved.
#' @param modification If an observation \code{xiim} is passed then \code{modification} will select the balancing method that \code{ccvc} uses to estimate the centred covariance.

#' @return If \code{boxes} is a list of numerical values then MVL is estimated for square boxes with side length given by \code{boxes}.
#'  The returned object is then an \code{fv} object containing estimates of MVL, box mass variance and box mass mean.
#'  If \code{boxes} is a list of owin objects then \code{mvlcc} returns a dataframe of with columns corresponding to estimates of MVL, box mass variance and box mass mean.
#'  Note if NA or NaN values in the \code{covariance} object are used then \code{mvlc} will return NA or NaN instead of an MVL value. 
#'  If the true covariance function and coverage probability of a RACS are passed to \code{mvlc} then the results will be the true MVL, box mass variance and box mass mean for the RACS.

#' @examples
#' xi <- heather$coarse
#' cencovar <- ccvcs(xi, obswin = Frame(xi), modifications = c("pickaH"))$pickaH
#' p <- area(xi) / area(Frame(xi))
#' sidelengths <- seq(0.3, 14, by = 0.2)
#' mvlccest <- mvlcc(sidelengths, cencovar, p))
#' # what is the MVL estimates for boxes that are discs?
#' discboxes <- lapply(sidelengths / 2, disc)
#' discmvls <- mvlcc(discboxes, cencovar, p)
#' # points(sidelengths, discmvls)
#' 
#' #direct to an image
#' xiim <- as.im(xi, na.replace = 0)
#' mvlccest <- mvlcc(sidelengths, xiim = xiim, modification = "pickaH")
#' 
#' @keywords spatial nonparametric 
mvlcc <- function(boxes, cencovar = NULL, p = NULL, xiim = NULL, modification = "pickaH"){
  if (!(is.null(cencovar) && is.null(p))){
    if (!is.null(xiim)){stop("xiim (an observation image) and cencovar or p were given. Either cencovar and p must be supplied or xiim supplied.")}
    lacv <- mvlcc.inputcovar(boxes, cencovar, p)
    unitname <- unitname(cencovar)
  } else if (!is.null(xiim)){
    p <- sum(xiim) / sum(is.finite(xiim$v))
    cencovar <- ccvcs(xiim, modifications = c(modification))[[1]]
    lacv <- mvlcc.inputcovar(boxes, cencovar, p)
    unitname <- unitname(xiim)
  } else {
    stop("Input requires specification of xiim or cencovar and p")
  }

  lacsdf <- as.data.frame(lacv)
  if (mode(boxes) %in% c("integer", "numeric")){
    lacsdf <- cbind(s = boxes, lacsdf)
    #recommended xlim:
    alim.min <- 1
    alim.max <- min(which(vapply(lacsdf[, "MVL"], is.na, FUN.VALUE = TRUE)), nrow(lacsdf))
    lacfv <- fv(lacsdf,
                argu = "s",
                valu = "MVL",
                fmla = ".y ~ s",
                alim = c(lacsdf[alim.min, "s"], lacsdf[alim.max, "s"]),
                ylab = quote(MVL[c]),
                unitname = unitname(xiim),
                labl = c("Box Width",
                         "MVL",
                         "Var(BoxMass)",
                         "Mean[BoxMass]"),
                desc = c("side lengths of boxes", 
                         "MVL derived from covariance",
                         "A covariance-based estimate of the variance in box mass",
                         "An estimate of mean box mass using coverage probability"
                )
    )
    fvnames(lacfv, a = ".") <- "MVL"
    return(lacfv)
  } else (return(lacsdf))
}

mvlcc.inputcovar <- function(boxes, cencovar, p){
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

  integrationresults <- mapply(innerprod.im, boxcov, list(cencovar), outsideA = 0, outsideB = NA, na.rm = FALSE, SIMPLIFY = FALSE) # the list around the cencovar is necessary to stop mapply unlisting the image itself

  lac <- unlist(integrationresults) / (p ^ 2 * boxarea ^ 2)
  return(list(
    MVL = lac,
    s2 = unlist(integrationresults),
    xbar = p * boxarea
  ))
}
