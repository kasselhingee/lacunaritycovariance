#' @title Pair-correlation based estimates of gliding box lacunarity
#' @export gblg
#'
#' @description Estimates the gliding box lacunarity (GBL) of a stationary RACS by estimating pair-correlation from a binary map (Hingee et al., 2017).
#'  It can also calculate the GBL of a RACS from a provided pair-correlation function. 

#' @references
#' Hingee K, Baddeley A, Caccetta P, Nair G (2019). Computation of lacunarity from covariance of spatial binary maps. \emph{Journal of Agricultural, Biological and Environmental Statistics}, 24, 264-288. DOI: 10.1007/s13253-019-00351-9.

#' @details
#' If we denote the estimated pair-correlation by \eqn{\hat{g}(v)}{g(v)} then the estimate of GBL is
#' \deqn{\frac{1}{|B|^2}\int \gamma_B(v)\hat{g}(v)dv, }{  \int gammaB(v) g(v) dv  /  (|B|^2),  }
#' where \eqn{B} is each of the sets (often called a box) specified by \code{boxes},
#' \eqn{\gamma_B}{gammaB} is the set covariance of \eqn{B},
#' \eqn{|B|} is the area of \eqn{B},
#' \eqn{p} is the coverage probability of a stationary RACS.
#' This can be used to compute the GBL from model parameters by passing \code{gblc} the 
#' covariance and coverage probability of the model.
#'
#' If the \code{xiim} argument to \code{gblg} is used then pair correlation is estimated from \code{xiim} using \code{\link{paircorr}} and the \code{pickaH} estimator.
#' 
#' The set covariance of \eqn{B} is computed empirically using \pkg{spatstat}'s \code{\link[spatstat]{setcov}} function, which converts \eqn{B} into a binary pixel mask using \code{\link[spatstat]{as.mask}} defaults. Computation speed can be increased by setting a small default number of pixels, \code{npixel}, in \pkg{spatstat}'s global options (accessed through \code{\link[spatstat]{spatstat.options}}), however fewer pixels also decreases the accuracy of the GBL computation.
#' 
#' The default integration method for this function uses [cubature::cubintegrate()] from the \pkg{cubature} package.
#' The 'harmonisesum' integration method is known to produce numerical artefacts (Section 6.2 of (Hingee et al., 2019))


#' @param boxes Either a list of side lengths for square boxes or a list of \code{owin} objects of any shape.
#' @param paircorr  A \code{im} object containing the pair-correlation function
#' @param xiim An observation of a stationary RACS as an \code{im} object. \code{xiim} must have values of either 1, 0 or NA; 1 denotes inside the RACS, 0 denotes outside, and NA denotes unobserved.
#' @param integrationMethod The integration method used by [innerprod.im()].

#' @return If \code{boxes} is a list of numerical values then GBL is estimated for square boxes with side length given by \code{boxes}.
#'  The returned object is then an \code{fv} object containing estimates of GBL.
#'  If \code{boxes} is a list of \code{owin} objects then \code{gblg} returns a dataframe of with columns corresponding to estimates of GBL.
#'  
#'  Note that if any values in the \code{paircorr} object needed for \code{gblg} are \code{NA} or \code{NaN} then \code{gblg} will return \code{NA} or \code{NaN}, respectively. 

#' @examples
#' xi <- as.im(heather$coarse, na.replace = 0, eps = 4 * heather$coarse$xstep)
#' sidelengths <- seq(0.3, 14, by = 3)
#'
#' # reduce resolution in setcov() for faster (less accurate) computation 
#' oldopt <- spatstat.options()
#' spatstat.options("npixel" = 2^4)
#' 
#' # compute GBL estimates from binary map
#' xiim <- as.im(xi, na.replace = 0)
#' gblgest <- gblg(sidelengths, xiim = xiim)
#'
#' spatstat.options(oldopt)
#' 
#' @keywords spatial nonparametric 
gblg <- function(boxes, paircorr = NULL, xiim = NULL, integrationMethod = "cubature"){
  if (!(is.null(paircorr))){
    if (!is.null(xiim)){stop("xiim (an observation image) and paircorr were given. paircorr and xiim cannot be simultaneously supplied.")}
    lacv <- gblg.inputpaircorr(boxes, paircorr)
    unitname <- unitname(paircorr)
  } else if (!is.null(xiim)){
    paircorr <- paircorr(xiim, estimators = "pickaH", drop = TRUE)
    lacv <- gblg.inputpaircorr(boxes, paircorr, integrationMethod = integrationMethod)
    unitname <- unitname(xiim)
  } else {
    stop("Input requires specification of xiim or paircorr")
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

gblg.inputpaircorr <- function(boxes, paircorr, integrationMethod = "cubature"){
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

  integrationresults <- mapply(innerprod.im, boxcov, list(paircorr - 1), outsideA = 0, outsideB = NA, na.replace = TRUE, method = integrationMethod, SIMPLIFY = FALSE) # the list around the paircorr is necessary to stop mapply unlisting the image itself

  GBLest <- 1 + unlist(integrationresults) / (boxarea ^ 2) 
  return(GBLest)
}


