#' @title Covariance-based calculations of mass variance lacunarity
#' @export lac
#' @export mvl
#' @export mvlc
#'
#' @description Estimates the mass variance lacunarity (MVL) of a stationary RACS from a bi-tonal image using the covariance method.
#'  It can also calculate the MVL of a RACS from a provided covariance (two-point probability) and coverage probability. 

#' @details
#' \code{mvl} and \code{lac} are aliases of mvlc.
#' If we denoted the estimated covariance by \eqn{\hat{C}(v)} and coverage probability \eqn{\hat{p}} then the estimate of MVL is
#' \deqn{\frac{1}{\hat{p}^2 |B|^2}\int \gamma_B(v)\hat{C}(v)dv -1 }

#' @param boxes Either a list of sidelengths for square boxes or a list of \code{owin} objects of any shape.
#' @param covariance  A \code{im} object containing the covariance function
#' @param p The coverage probability. Typically estimated by the fraction of the observation window covered by the set of interest.
#' @param xiim An observation of a stationary RACS in \code{im} format. \code{xiim} must have values of either 1, 0 or NA; 1 denotes inside the RACS, 0 denotes outside, and NA denotes unobserved.

#' @return Either an \code{fv} object containing the MVL values and box side lengths or if \code{boxes} is a list of owin objects then \code{mvlc} returns a list of corresponding MVL values.
#'  Note if NA or NaN values in the \code{covariance} object are used then \code{mvlc} will return NA or NaN instead of an MVL value. 

#' @examples
#' xi <- heather$coarse
#' covar <- racscovariance(xi, inclraw = FALSE)
#' p <- area(xi) / area(Frame(xi))
#' sidelengths <- seq(0.3, 14, by = 0.2)
#' plot(lac(sidelengths, covar, p))
#' # what is the MVL estimates for boxes that are discs?
#' discboxes <- lapply(sidelengths / 2, disc)
#' discmvls <- mvlc(discboxes, covar, p)
#' points(sidelengths, discmvls)
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
    covar <- racscovariance(xiim, obswin = w)
    lacv <- mvlc.inputcovar(boxes, covar, p)
    unitname <- unitname(xiim)
  } else {
    stop("Input requires specification of xiim or covariance and p")
  }

  if (mode(boxes) %in% c("integer", "numeric")){
    lacfv <- fv(data.frame(s = boxes, MVL = lacv),
                argu = "s",
                valu = "MVL",
                ylab = expression(MVL),
                unitname = unitname,
                labl = c("Box Side Length", "MVL"),
                desc = c("Side length of boxes", "MVL derived from covariance")
               )
    return(lacfv)
  } else (return(lacv))
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

  integrationresults <- mapply(innerprod.im, boxcov, list(covariance), na.rm = FALSE, SIMPLIFY = FALSE) # the list around the covariance is necessary to stop mapply unlisting the image itself

  lac <- unlist(integrationresults) / (p ^ 2 * boxarea ^ 2) - 1
  return(lac)
}


innerprod.im <- function(A, B, na.rm = FALSE){
  integrationregion <- intersect.owin(Frame(A), Frame(B))
  #got to do the harmonisation manually so that NA values that the subsetting operation doesn't introduce NA values
  harmgrid <- as.mask(integrationregion,
             eps = c(min(A$xstep, B$xstep), min(A$ystep, B$ystep)))
  A2 <- as.im(A, xy = harmgrid)
  B2 <- as.im(B, xy = harmgrid)
  prdimg <- eval.im(A2 * B2, harmonize = FALSE)
  return(sum(prdimg[, ], na.rm = na.rm) * prdimg$xstep * prdimg$ystep)
}
#tests of innerprod.im:
#innerprod.im(as.im(square(1)),as.im(square(1),value=1))
#natest:
#imna <- as.im(square(1.01),value=1)
#imna[as.ppp(c(0.5,0.5),W=Frame(imna))] <- NA
#innerprod.im(as.im(square(1)),imna, na.rm=FALSE)
#innerprod.im(as.im(square(1)),imna, na.rm=TRUE)

#innerprod.im(as.im(square(1)),as.im(square(0.25),value=1))

#should be close to 0 (orthogonal):
#innerprod.im(as.im(function(x,y) {sin(x)},W=square(7*pi),eps=0.01),as.im(function(x,y) {sin(2*x)},W=square(2*pi),eps=0.01))
#should be very non-zero
#innerprod.im(as.im(function(x,y) {sin(x)},W=square(7*pi),eps=0.01),as.im(function(x,y) {sin(x)},W=square(2*pi),eps=0.01))
#it should be (and is) equal to this: sum(as.im(function(x,y){sin(x)*sin(x)},W=square(2*pi),eps=0.01))*0.01*0.01


#' @rdname mvlc 
mvl <- mvlc

#' @rdname mvlc
lac <- mvlc
