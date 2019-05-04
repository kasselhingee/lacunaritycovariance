#' @title Gliding box lacunarity estimator using plug-in moment covariance estimator
#' @export gblc
#'
#' @description 
#' Can be used to estimate the gliding box lacunarity (GBL) of a stationary RACS from a binary map
#'  using the plug-in moment covariance covariance estimator (Hingee et al., 2017).
#'  It can also calculate the GBL of a RACS from a given covariance function and coverage probability. 

#' @references
Hingee K, Baddeley A, Caccetta P, Nair G (2019). Computation of lacunarity from covariance of spatial binary maps. \emph{Journal of Agricultural, Biological and Environmental Statistics}, 24, 264-288. DOI: 10.1007/s13253-019-00351-9.

#' @details
#' Computes a numerical approximation of 
#' \deqn{\int \gamma_B(v) C(v) dv / (p^2 |B|^2).}{\int gammaB(v) C(v) dv / (p^2 |B|^2),}
#' where \eqn{B} is a given set (often called a box),
#' \eqn{\gamma_B}{gammaB} is the set covariance of \eqn{B},
#' \eqn{|B|} is the area of \eqn{B},
#' \eqn{p} is the coverage probability of a stationary RACS, and
#' \eqn{C(v)} is the covariance of a stationary RACS.
#' This can be used to compute the GBL from model parameters by passing \code{gblc} the 
#' covariance and coverage probability of the model.
#' 
#' If a binary map is supplied then \eqn{p} and \eqn{C(v)} are estimated using
#'  the usual coverage probability estimator and the plug-in moment covariance estimator, respectively 
#'  (see \code{\link{coverageprob}} and \code{\link{plugincvc}}).

#' @param boxes Either a list of sidelengths for square boxes or a list of \code{owin} objects of any shape.
#' @param covariance  A \code{im} object containing the covariance function
#' @param p The coverage probability. Typically estimated by the fraction of the observation window covered by the set of interest.
#' @param xiim An observation of a stationary RACS in \code{im} format. \code{xiim} must have values of either 1, 0 or NA; 1 denotes inside the RACS, 0 denotes outside, and NA denotes unobserved.

#' @return If \code{boxes} is a list of numerical values then GBL is estimated 
#' for square boxes with side length given by \code{boxes}.
#'  The returned object is then an \code{fv} object containing estimates of GBL,
#'   box mass variance and box mass mean.
#'  If \code{boxes} is a list of owin objects then \code{gblc} returns a 
#'  dataframe of with columns corresponding to estimates of GBL, box mass variance and box mass mean.
#' 
#'   Note if NA or NaN values in the \code{covariance} object are used then \code{gblc} will return NA or NaN. 

#' @examples
#' xi <- heather$coarse
#' covar <- plugincvc(xi, Frame(xi))
#' p <- area(xi) / area(Frame(xi))
#' if(interactive()){
#' sidelengths <- seq(0.3, 14, by = 0.2)
#' } else {
#' sidelengths <- seq(0.3, 14, by = 1)
#' }
#' gblest <- gblc(sidelengths, covar, p)
#' # plot(gblest)
#' # what is the GBL estimates for boxes that are discs?
#' discboxes <- lapply(sidelengths / 2, disc)
#' discgbls <- gblc(discboxes, covar, p)
#' # points(sidelengths, discgbls$GBL)
#' 
#' @keywords spatial nonparametric 
gblc <- function(boxes, covariance = NULL, p = NULL, xiim = NULL){
  if (!(is.null(covariance) && is.null(p))){
    if (!is.null(xiim)){stop("xiim (an observation image) and covariance or p were given. Either covariance and p must be supplied or xiim supplied.")}
    lacv <- gblc.inputcovar(boxes, covariance, p)
    unitname <- unitname(covariance)
  } else if (!is.null(xiim)){
    p <- sum(xiim) / sum(is.finite(xiim$v))
    w <- as.owin(xiim) #w is observation window - only the non NA values end up in window
    xiim[is.na(xiim$v)] <- 0
    covar <- plugincvc(xiim, obswin = w)
    lacv <- gblc.inputcovar(boxes, covar, p)
    unitname <- unitname(xiim)
  } else {
    stop("Input requires specification of xiim or covariance and p")
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
                unitname = unitname,
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

gblc.inputcovar <- function(boxes, covariance, p){
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

  GBLest <- unlist(integrationresults) / (p ^ 2 * boxarea ^ 2) 
  return(list(
    GBL = GBLest,
    s2 = unlist(integrationresults) -  (p ^ 2 * boxarea ^ 2),
    xbar = p * boxarea
  ))
}

