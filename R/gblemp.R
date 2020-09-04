#' @title Empirical Gliding Box Lacunarity
#' @export gblemp  gbltrad
#'
#' @description Calculates empirical gliding box lacunarity of a binary map, which was proposed by Allain and Cloitre (1991). 
#' @details Calculates empirical gliding box lacunarity (Allain and Cloitre, 1991) for a given range of square box sizes,
#' \deqn{1 + Var(area(B . xi)) / E[area(B . xi)]^2, }
#' where \eqn{B} is a box that has a random location in the observation window and \eqn{area(B . xi)} is the (random) area of the foreground in \eqn{B}. 
#' This is an estimate of the gliding box lacunarity of a RACS (Hingee et al., 2017).
#' 
#' The algorithm uses the pixel locations in \code{xiim} as an array of box centre locations to compute
#'  the mean and variance of the area in a random box of a given size.
#' Locations where the box is not completely within the observation window are ignored.
#'  
#' @section WARNING: 
#' The box side lengths are rounded such that they are an odd number of pixels across.
#' \code{gblemp} uses the \code{\link[RcppRoll]{roll_sum}} function in \pkg{RcppRoll} to operate, so \pkg{RcppRoll} must be installed to run \code{gblemp}.



#' @references 
#' Allain, C. and Cloitre, M. (1991) Characterizing the lacunarity of random and deterministic fractal sets. \emph{Physical Review A}, 44, 3552-3558.
#' 
#' Hingee K, Baddeley A, Caccetta P, Nair G (2019). Computation of lacunarity from covariance of spatial binary maps. \emph{Journal of Agricultural, Biological and Environmental Statistics}, 24, 264-288. DOI: 10.1007/s13253-019-00351-9.
#'
#' @return An \code{fv} object containing empirical GBL, variance of the area in the box and mean of the area in the box. 
#'  The box widths (labelled \code{s}) are always odd multiples of the pixel width.
#'  
#' @param xiim An image of pixels valued either \code{0}, \code{1} or \code{NA}. \code{NA} valued pixels are assumed to be outside the observation window.
#' @param boxwidths A list of suggested box widths in the same units as \code{xiim}. 
#' Note the actual box widths used by \code{gblemp} will be the closest multiple of an odd number of pixel widths.
#' @param obswin Optional observation window. The observation window used for the estimator will be the intersection of \code{obswin} and the pixels that are not \code{NA} in \code{xiim}.
#' @examples
#' xiim <- as.im(heather$coarse, na.replace = 0)
#' boxwidths <- seq(0.2, 14, by = 0.2) #in units of xiim
#' gblest <- gblemp(boxwidths, xiim)
#'
#' @keywords spatial nonparametric 
gblemp <- function(boxwidths, xiim, obswin = Frame(xiim)){
  if (!is.im(xiim)){stop("input xiim must be of class im")}
  if (abs(xiim$xstep - xiim$ystep) > 1E-2 * xiim$xstep){stop("image pixels must be square")}
  isbinarymap(xiim, requiretrue = TRUE)
#convert boxwidths to odd pixel amounts, taking into account that want a distance to edge
  spix <- 1 + round( (boxwidths - xiim$xstep) / (2 * xiim$xstep)) * 2
  spix <- unique(spix)
  rpix <- (spix - 1) / 2
  sidel <- spix * xiim$xstep

#compute observation mask
  obsvd <- xiim
  obsvd[is.finite(xiim$v)] <- TRUE
  if (class(obswin) == "im"){obsvd <- eval.im(obswin * obsvd)}
  obsvd <- as.owin(obsvd) #owin format may not be needed anymore
  if (class(obswin) == "owin"){obsvd <- intersect.owin(obsvd, obswin)}

  if (requireNamespace("RcppRoll") != TRUE){
     stop("RcppRoll package must be installed to calculate empirical gliding box lacunarity")
  }

xiim[(complement.owin(intersect.owin(obswin, Frame(xiim)), frame = Frame(xiim)))] <- NA  #make sure the pixels outside obswin are set to NA so that reduce sampling happens naturally ##NOTE: this a time consuming operation that may never be needed
lacs <- mapply(gblemp_intern.rcpproll, sidep = 2 * rpix + 1, MoreArgs = list(xiim = xiim, obswin = obsvd), SIMPLIFY = FALSE)

  #converting results in fv objects
    valsdf <- matrix(unlist(lacs), ncol = length(lacs[[1]]), byrow = TRUE)
    colnames(valsdf) <- names(lacs[[1]])
    lacsdf <- cbind(data.frame(s = sidel), valsdf)
    #recommended xlim:
    alim.min <- 1
    alim.max <- min(which(vapply(lacsdf[, "GBL"], is.na, FUN.VALUE = TRUE)), nrow(lacsdf))
    lacfv <- fv(lacsdf,
           argu = "s",
           valu = "GBL",
           fmla = ".y ~ s",
           alim = c(lacsdf[alim.min, "s"], lacsdf[alim.max, "s"]),
           ylab = expression(GBL[gb]),
           unitname = unitname(xiim),
           labl = c("Box Width",
                    "GBL",
                    "Var(BoxMass)",
                    "Avg[BoxMass]"),
           desc = c("side lengths of boxes", 
                    "Empirical GBL",
                    "Variance of foreground area in the randomly placed box in empirical GBL",
                    "Mean foreground area in the randomly placed box in empirical GBL"
                    ),
           fname = "GBL"
           )
    fvnames(lacfv, a = ".") <- "GBL"
  return(lacfv)
}


##########################
##The following function calculates lacunarity for a box with side lengths 2*bX+1 and 2*bY+1 (in pixels). The RS version is automatically calculated by ignoring those boxes that have sums that includa NA values. 
#eg gblemp_intern.rcpproll(xiim,5,5,5*0.8)
#the obswin is only for the raw version and must be an owin object. 
#uses rcpproll
gblemp_intern.rcpproll <- function(xiim, sidep, obswin = Frame(xiim)){
  mat <- as.matrix(xiim)
  if ( (sidep > nrow(mat)) | (sidep > ncol(mat))){
    gblemp.rs <- NA
    sampmean.rs <- NA
    samp2ndmom.rs <- NA
  }
  else {
    movline.overrows <- RcppRoll::roll_sum(mat, sidep)
    movline.overrowthencols <- RcppRoll::roll_sum(t(movline.overrows), sidep) * xiim$xstep * xiim$ystep
    sampmean.rs <- mean(movline.overrowthencols, na.rm = TRUE) #sample mean
    samp2ndmom.rs <- mean(movline.overrowthencols ^ 2, na.rm = TRUE) #biased sample second moment
    gblemp.rs <- samp2ndmom.rs / (sampmean.rs ^ 2) 
  }

    return(list(
      GBL = gblemp.rs,
      s2 = samp2ndmom.rs - sampmean.rs^2,
      xbar = sampmean.rs
      ))
}

#' @describeIn gblemp An alias of \code{gblemp} used in past versions of this package. This alias may be removed in future versions.
gbltrad <- function(boxwidths, xiim, obswin = Frame(xiim)){
  warning("'gbltrad' function name has been changed to 'gblemp'. Please use 'gblemp' instead of 'gbltrad'.")
  gblemp(boxwidths, xiim, obswin)
}

