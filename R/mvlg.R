#' @title Pair-correlation based calculations of mass variance lacunarity
#' @export mvlg
#'
#' @description Estimates the mass variance lacunarity (MVL) of a stationary RACS from a bi-tonal image by estimating pair-correlation.
#'  It can also calculate the MVL of a RACS from a provided pair-correlation function. 

#' @details
#' If we denoted the estimated pair-correlation by \eqn{\hat{g}(v)} then the estimate of MVL is
#' \deqn{\frac{1}{|B|^2}\int \gamma_B(v)\hat{g}(v)dv - 1 }

#' @param boxes Either a list of sidelengths for square boxes or a list of \code{owin} objects of any shape.
#' @param paircor  A \code{im} object containing the pair-correlation function
#' @param xiim An observation of a stationary RACS in \code{im} format. \code{xiim} must have values of either 1, 0 or NA; 1 denotes inside the RACS, 0 denotes outside, and NA denotes unobserved.

#' @return Either an \code{fv} object containing the MVL values and box side lengths or if \code{boxes} is a list of owin objects then \code{mvlg} returns a list of corresponding MVL values.
#'  Note if NA or NaN values in the \code{paircor} object are used then \code{mvlg} will return NA or NaN instead of an MVL value. 

#' @examples
#' xi <- heather$coarse
#' paircor <- pclns(as.im(xi, na.replace = 0), modifications = "pickahajek")[[1]]
#' sidelengths <- seq(0.3, 14, by = 0.2)
#' plot(mvlg(sidelengths, paircor))
#' # what is the MVL estimates for boxes that are discs?
#' discboxes <- lapply(sidelengths / 2, disc)
#' discmvls <- mvlg(discboxes, paircor)
#' points(sidelengths, discmvls)
#' 
#' @keywords spatial nonparametric 
mvlg <- function(boxes, paircor = NULL, xiim = NULL){
  if (!(is.null(paircor))){
    if (!is.null(xiim)){stop("xiim (an observation image) and paircor were given. paircor and xiim cannot be simultaneously supplied.")}
    lacv <- mvlg.inputpaircor(boxes, paircor)
    unitname <- unitname(paircor)
  } else if (!is.null(xiim)){
    paircor <- pclns(xiim, modifications = "pickahajek")[[1]]
    lacv <- mvlg.inputpaircor(boxes, paircor)
    unitname <- unitname(xiim)
  } else {
    stop("Input requires specification of xiim or paircor and p")
  }

  if (mode(boxes) %in% c("integer", "numeric")){
    lacfv <- fv(data.frame(s = boxes, MVL = lacv),
                argu = "s",
                valu = "MVL",
                ylab = expression(MVL),
                unitname = unitname,
                labl = c("Box Side Length", "MVL"),
                desc = c("Side length of boxes", "MVL derived from pair-correlation")
               )
    return(lacfv)
  } else (return(lacv))
}

mvlg.inputpaircor <- function(boxes, paircor){
  stopifnot(is.im(paircor))
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

  integrationresults <- mapply(innerprod.im, boxcov, list(paircor), na.rm = FALSE, SIMPLIFY = FALSE) # the list around the paircor is necessary to stop mapply unlisting the image itself

  lac <- unlist(integrationresults) / (boxarea ^ 2) - 1
  return(lac)
}


