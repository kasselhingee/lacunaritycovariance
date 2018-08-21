#' @title Pair-correlation based estimates of mass variance lacunarity
#' @export mvlg
#'
#' @description Estimates the mass variance lacunarity (MVL) of a stationary RACS by estimating pair-correlation from a binary map.
#'  It can also calculate the MVL of a RACS from a provided pair-correlation function. 

#' @details
#' If we denoted the estimated pair-correlation by \eqn{\hat{g}(v)}{g(v)} then the estimate of MVL is
#' \deqn{\frac{1}{|B|^2}\int \gamma_B(v)\hat{g}(v)dv - 1. }{  \int gammaB(v) g(v) dv  /  (|B|^2)  - 1.  }

#' @param boxes Either a list of sidelengths for square boxes or a list of \code{owin} objects of any shape.
#' @param paircorr  A \code{im} object containing the pair-correlation function
#' @param xiim An observation of a stationary RACS in \code{im} format. \code{xiim} must have values of either 1, 0 or NA; 1 denotes inside the RACS, 0 denotes outside, and NA denotes unobserved.

#' @return If \code{boxes} is a list of numerical values then MVL is estimated for square boxes with side length given by \code{boxes}.
#'  The returned object is then an \code{fv} object containing estimates of MVL.
#'  If \code{boxes} is a list of owin objects then \code{mvlg} returns a dataframe of with columns corresponding to estimates of MVL.
#'  
#'  Note if NA or NaN values in the \code{paircorr} object are used then \code{mvlg} will return NA or NaN instead of an MVL value. 

#' @examples
#' xi <- heather$coarse
#' pcln <- paircorr(as.im(xi, na.replace = 0), modifications = "pickaH", drop = TRUE)
#' sidelengths <- seq(0.3, 14, by = 0.2)
#' mvlgest <- mvlg(sidelengths, pcln)
#' # what is the MVL estimates for boxes that are discs?
#' discboxes <- lapply(sidelengths / 2, disc)
#' discmvls <- mvlg(discboxes, pcln)
#' # points(sidelengths, discmvls)
#' 
#' @keywords spatial nonparametric 
mvlg <- function(boxes, paircorr = NULL, xiim = NULL){
  if (!(is.null(paircorr))){
    if (!is.null(xiim)){stop("xiim (an observation image) and paircorr were given. paircorr and xiim cannot be simultaneously supplied.")}
    lacv <- mvlg.inputpaircorr(boxes, paircorr)
    unitname <- unitname(paircorr)
  } else if (!is.null(xiim)){
    paircorr <- paircorr(xiim, modifications = "pickaH", drop = TRUE)
    lacv <- mvlg.inputpaircorr(boxes, paircorr)
    unitname <- unitname(xiim)
  } else {
    stop("Input requires specification of xiim or paircorr and p")
  }

  if (mode(boxes) %in% c("integer", "numeric")){
    lacfv <- fv(data.frame(s = boxes, MVL = lacv),
                argu = "s",
                valu = "MVL",
                ylab = expression(MVL),
                unitname = unitname,
                labl = c("Box Side Length", "MVL"),
                desc = c("Side length of boxes", "MVL derived from pair-correlation"),
                fname = "MVL"
               )
    return(lacfv)
  } else (return(lacv))
}

mvlg.inputpaircorr <- function(boxes, paircorr){
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

  lac <- unlist(integrationresults) / (boxarea ^ 2) - 1
  return(lac)
}


