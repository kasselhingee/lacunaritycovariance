#' @title Pair-correlation based estimates of gliding box lacunarity
#' @export gblg
#'
#' @description Estimates the gliding box lacunarity (GBL) of a stationary RACS by estimating pair-correlation from a binary map (Hingee et al., 2017).
#'  It can also calculate the GBL of a RACS from a provided pair-correlation function. 

#' @references
#' Hingee K, Baddeley A, Caccetta P, Nair G (2019). Computation of lacunarity from covariance of spatial binary maps. \emph{Journal of Agricultural, Biological and Environmental Statistics}, 24, 264-288. DOI: 10.1007/s13253-019-00351-9.

#' @details
#' If we denote the estimated pair-correlation by \eqn{\hat{g}(v)}{g(v)} then the estimate of GBL is
#' \deqn{\frac{1}{|B|^2}\int \gamma_B(v)\hat{g}(v)dv. }{  \int gammaB(v) g(v) dv  /  (|B|^2) .  }

#' @param boxes Either a list of side lengths for square boxes or a list of \code{owin} objects of any shape.
#' @param paircorr  A \code{im} object containing the pair-correlation function
#' @param xiim An observation of a stationary RACS as an \code{im} object. \code{xiim} must have values of either 1, 0 or NA; 1 denotes inside the RACS, 0 denotes outside, and NA denotes unobserved.

#' @return If \code{boxes} is a list of numerical values then GBL is estimated for square boxes with side length given by \code{boxes}.
#'  The returned object is then an \code{fv} object containing estimates of GBL.
#'  If \code{boxes} is a list of \code{owin} objects then \code{gblg} returns a dataframe of with columns corresponding to estimates of GBL.
#'  
#'  Note if value in the \code{paircorr} object that are needed for \code{gblg} are \code{NA} or \code{NaN} then \code{gblg} will return \code{NA} or \code{NaN}. 

#' @examples
#' xi <- heather$coarse
#' pcln <- paircorr(as.im(xi, na.replace = 0), estimators = "pickaH", drop = TRUE)
#' if (interactive()) {
#' sidelengths <- seq(0.3, 14, by = 0.2)
#' } else {
#' sidelengths <- seq(0.3, 14, by = 1)
#' }
#' gblgest <- gblg(sidelengths, pcln)
#' # what is the GBL estimates for boxes that are discs?
#' discboxes <- lapply(sidelengths / 2, disc)
#' discgbls <- gblg(discboxes, pcln)
#' # points(sidelengths, discgbls)
#' 
#' @keywords spatial nonparametric 
gblg <- function(boxes, paircorr = NULL, xiim = NULL){
  if (!(is.null(paircorr))){
    if (!is.null(xiim)){stop("xiim (an observation image) and paircorr were given. paircorr and xiim cannot be simultaneously supplied.")}
    lacv <- gblg.inputpaircorr(boxes, paircorr)
    unitname <- unitname(paircorr)
  } else if (!is.null(xiim)){
    paircorr <- paircorr(xiim, estimators = "pickaH", drop = TRUE)
    lacv <- gblg.inputpaircorr(boxes, paircorr)
    unitname <- unitname(xiim)
  } else {
    stop("Input requires specification of xiim or paircorr and p")
  }

  if (mode(boxes) %in% c("integer", "numeric")){
    lacfv <- fv(data.frame(s = boxes, GBL = lacv),
                argu = "s",
                valu = "GBL",
                ylab = expression(GBL),
                unitname = unitname,
                labl = c("Box Side Length", "GBL"),
                desc = c("Side length of boxes", "GBL derived from pair-correlation"),
                fname = "GBL"
               )
    return(lacfv)
  } else (return(lacv))
}

gblg.inputpaircorr <- function(boxes, paircorr){
  stopifnot(is.im(paircorr))
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

  integrationresults <- mapply(innerprod.im, boxcov, list(paircorr), outsideA = 0, outsideB = NA, na.rm = FALSE, SIMPLIFY = FALSE) # the list around the paircorr is necessary to stop mapply unlisting the image itself

  GBLest <- unlist(integrationresults) / (boxarea ^ 2) 
  return(GBLest)
}


