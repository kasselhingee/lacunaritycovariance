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

#' @return Either an \code{fv} object containing the MVL values and box side lengths or if \code{boxes} is a list of owin objects then \code{mvlc} returns a list of corresponding MVL values.
#'  Note if NA or NaN values in the \code{cencovar} object are used then \code{mvlc} will return NA or NaN instead of an MVL value. 

#' @examples
#' xi <- heather$coarse
#' cencovar <- ccvcs(xi, obswin = Frame(xi), modifications = c("pickaH"))$pickaH
#' p <- area(xi) / area(Frame(xi))
#' sidelengths <- seq(0.3, 14, by = 0.2)
#' # plot(mvlcc(sidelengths, cencovar, p))
#' # what is the MVL estimates for boxes that are discs?
#' discboxes <- lapply(sidelengths / 2, disc)
#' discmvls <- mvlcc(discboxes, cencovar, p)
#' # points(sidelengths, discmvls)
#' 
#' #direct to an image
#' xiim <- as.im(xi, na.replace = 0)
#' # plot(mvlcc(sidelengths, xiim = xiim, modification = "pickaH"))
#' # plot(add = TRUE, mvlc(sidelengths, xiim = xiim), lty = "dashed", col = "red")
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

  if (mode(boxes) %in% c("integer", "numeric")){
    lacfv <- fv(data.frame(s = boxes, MVL = lacv),
                argu = "s",
                valu = "MVL",
                ylab = expression(MVL),
                unitname = unitname,
                labl = c("Box Side Length", "MVL"),
                desc = c("Side length of boxes", "MVL derived from centred covariance")
               )
    return(lacfv)
  } else (return(lacv))
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
  return(lac)
}
