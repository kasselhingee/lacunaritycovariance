#' @title Covariance-based calculations of mass variance lacunarity
#' @export mvlc
#'
#' @description Estimates the mass variance lacunarity (MVL) of a stationary RACS from a bi-tonal image using the covariance method.
#'  It can also calculate the MVL of a RACS from a provided covariance (two-point probability) and coverage probability. 

#' @details
#' If we denoted the estimated covariance by \eqn{\hat{C}(v)} and coverage probability \eqn{\hat{p}} then the estimate of MVL is
#' \deqn{\frac{1}{\hat{p}^2 |B|^2}\int \gamma_B(v)\hat{C}(v)dv -1 }

#' @param boxes Either a list of sidelengths for square boxes or a list of \code{owin} objects of any shape.
#' @param covariance  A \code{im} object containing the covariance function
#' @param p The coverage probability. Typically estimated by the fraction of the observation window covered by the set of interest.
#' @param xiim An observation of a stationary RACS in \code{im} format. \code{xiim} must have values of either 1, 0 or NA; 1 denotes inside the RACS, 0 denotes outside, and NA denotes unobserved.

#' @return If \code{boxes} is a list of numerical values then MVL is estimated for square boxes with side length given by \code{boxes}.
#'  The returned object is then an \code{fv} object containing estimates of MVL, box mass variance and box mass mean.
#'  If \code{boxes} is a list of owin objects then \code{mvlc} returns a dataframe of with columns corresponding to estimates of MVL, box mass variance and box mass mean.
#'  Note if NA or NaN values in the \code{covariance} object are used then \code{mvlc} will return NA or NaN instead of an MVL value. 
#'  If the true covariance function and coverage probability of a RACS are passed to \code{mvlc} then the results will be the true MVL, box mass variance and box mass mean for the RACS.

#' @examples
#' xi <- heather$coarse
#' covar <- tradcovarest(xi, inclraw = FALSE)
#' p <- area(xi) / area(Frame(xi))
#' sidelengths <- seq(0.3, 14, by = 0.2)
#' mvlest <- mvlc(sidelengths, covar, p)
#' # plot(mvlest)
#' # what is the MVL estimates for boxes that are discs?
#' discboxes <- lapply(sidelengths / 2, disc)
#' discmvls <- mvlc(discboxes, covar, p)
#' # points(sidelengths, discmvls$MVL)
#' 
#' @keywords spatial nonparametric 
mvlc <- function(boxes, covariance = NULL, p = NULL, xiim = NULL){
  if (!(is.null(covariance) && is.null(p))){
    if (!is.null(xiim)){stop("xiim (an observation image) and covariance or p were given. Either covariance and p must be supplied or xiim supplied.")}
    lacv <- mvlc.inputcovar(boxes, covariance, p)
    unitname <- unitname(covariance)
  } else if (!is.null(xiim)){
    p <- sum(xiim) / sum(is.finite(xiim$v))
    w <- as.owin(xiim) #w is observation window - only the non NA values end up in window
    xiim[is.na(xiim$v)] <- 0
    covar <- tradcovarest(xiim, obswin = w)
    lacv <- mvlc.inputcovar(boxes, covar, p)
    unitname <- unitname(xiim)
  } else {
    stop("Input requires specification of xiim or covariance and p")
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
                unitname = unitname,
                labl = c("Box Width",
                         "MVL",
                         "Var(BoxMass)",
                         "Mean[BoxMass]"),
                desc = c("side lengths of boxes", 
                         "MVL derived from covariance",
                         "A covariance-based estimate of the variance in box mass",
                         "An estimate of mean box mass using coverage probability"
                ),
                fname = "MVL"
    )
    fvnames(lacfv, a = ".") <- "MVL"
    return(lacfv)
  } else (return(lacsdf))
}

mvlc.inputcovar <- function(boxes, covariance, p){
  stopifnot(is.im(covariance))
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

  integrationresults <- mapply(innerprod.im, boxcov, list(covariance), outsideA = 0, outsideB = NA, na.rm = FALSE, SIMPLIFY = FALSE) # the list around the covariance is necessary to stop mapply unlisting the image itself

  lac <- unlist(integrationresults) / (p ^ 2 * boxarea ^ 2) - 1
  return(list(
    MVL = lac,
    s2 = unlist(integrationresults) -  (p ^ 2 * boxarea ^ 2),
    xbar = p * boxarea
  ))
}

