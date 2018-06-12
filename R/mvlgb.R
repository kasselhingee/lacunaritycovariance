#' @title Gliding Box estimation of mass variance lacunarity
#' @export mvlgb lacgb
#' @importFrom utils installed.packages
#'
#' @description Calculates the gliding box estimate [1] of mass variance lacunarity from a bi-tonal image.
#' @details Calculates the gliding box estimate [1] of mass variance lacunarity for a given range of square box sizes. 
#' The algorithm uses the pixel locations in \code{xiim} as an array of box centre locations to compute
#'  the mean and variance of the area in a random box of a given size.
#' Locations where the box is not completely within the observation window are ignored.
#'  
#' WARNING: This function needs the \code{roll_sum} function in \code{RcppRoll} to operate.
#' 
#' Note: The side lengths are rounded such that they are an odd number of pixels across.
#' 
#' @references [1] Allain, C. and Cloitre, M. (1991) Characterizing the lacunarity of random and deterministic fractal sets. Physical Review A, 44, 3552-3558.
#'
#' @return An \code{fv} object with a column called \code{MVL} for the usual gliding box estimate described in [1]. 
#'  The side lengths (labelled \code{s}) are always odd multiples of the pixel width.
#'  
#' Another column, called \code{raw}, is included if \code{inclraw=TRUE}. 
#' This column represents a version of the gliding box algorithm that pretends that the items/cover of interest are fully observed
#'  (i.e. that the observation window is the entire plane).
#' @param xiim An image of pixels valued either \code{0}, \code{1} or \code{NA}. \code{NA} valued pixels are assumed to be outside the observation window.
#' @param sidelengths A list of suggested box side lengths in the same units as \code{xiim}. Note the actual side lengths used will be the closest multiple of an odd number of pixel widths.
#' @param inclraw If TRUE the function will also return a gliding box estimate that ignores edge effects. (Default is FALSE)
#' @param obswin Optional observation window. The observation window used for the estimator will be the intersection of \code{obswin} and the pixels that are not \code{NA} in \code{xiim}.
#' @examples
#' xiim <- as.im(heather$coarse, na.replace = 0)
#' sidelengths <- seq(0.2, 14, by = 0.2) #in units of xiim
#' lac <- mvlgb(sidelengths, xiim, inclraw = FALSE)
#' # plot(lac)
#'
#' @keywords spatial nonparametric 
mvlgb <- function(sidelengths, xiim, inclraw = FALSE, obswin = Frame(xiim)){
  if (!is.im(xiim)){stop("input xiim must be of class im")}
  if (abs(xiim$xstep - xiim$ystep) > 1E-2 * xiim$xstep){stop("image pixels must be square")}
#convert sidelengths to odd pixel amounts, taking into account that want a distance to edge
  spix <- 1 + round( (sidelengths - xiim$xstep) / (2 * xiim$xstep)) * 2
  spix <- unique(spix)
  rpix <- (spix - 1) / 2
  sidel <- spix * xiim$xstep

#compute observation mask
  obsvd <- xiim
  obsvd[is.finite(xiim$v)] <- TRUE
  if (class(obswin) == "im"){obsvd <- eval.im(obswin * obsvd)}
  obsvd <- as.owin(obsvd) #owin format may not be needed anymore
  if (class(obswin) == "owin"){obsvd <- intersect.owin(obsvd, obswin)}

  if (!("RcppRoll" %in% installed.packages()[, 1])){
     stop("RcppRoll must be installed to calculate gliding box lacunarity")
  }

xiim[(complement.owin(intersect.owin(obswin, Frame(xiim)), frame = Frame(xiim)))] <- NA  #make sure the pixels outside obswin are set to NA so that reduce sampling happens naturally ##NOTE: this a time consuming operation that may never be needed
lacs <- mapply(mvlgb_intern.rcpproll, sidep = 2 * rpix + 1, MoreArgs = list(xiim = xiim, obswin = obsvd), SIMPLIFY = FALSE, inclraw)

  #converting results in fv objects
  
  #extra terms if raw is exported
  if (inclraw){
    nobord <- unlist(lapply(lacs, `[[`, 1) )
    MVL <- unlist(lapply(lacs, `[[`, 2) )
    lacsdf <- data.frame(s = sidel, raw = nobord, MVL = MVL)
    lacfv <- fv(lacsdf, argu = "s", valu = "MVL",
           ylab = expression(MVL[gb]),
           unitname = unitname(xiim),
           labl = c("side length", "raw", "MVL"),
           desc = c("side lengths of boxes", "Gliding Box Lacunarity ignoring edge effects", "Gliding Box Lacunarity that only uses boxes entirely within the observation")
           )
  }
  else {
    valsdf <- matrix(unlist(lacs), ncol = length(lacs[[1]]), byrow = TRUE)
    colnames(valsdf) <- names(lacs[[1]])
    lacsdf <- cbind(data.frame(s = sidel), valsdf)
    #recommended xlim:
    alim.min <- 1
    alim.max <- min(which(vapply(lacsdf[, "MVL"], is.na, FUN.VALUE = TRUE)), nrow(lacsdf))
    lacfv <- fv(lacsdf,
           argu = "s",
           valu = "MVL",
           fmla = ".y ~ s",
           alim = c(lacsdf[alim.min, "s"], lacsdf[alim.max, "s"]),
           ylab = expression(MVL[gb]),
           unitname = unitname(xiim),
           labl = c("Box Width",
                    "MVL",
                    "Var(BoxMass)",
                    "Avg[BoxMass]"),
           desc = c("side lengths of boxes", 
                    "Gliding box MVL estimate that only uses boxes entirely within the observation",
                    "Variance of box mass - used in the gliding box estimate",
                    "Average box mass - used in the gliding box estimate"
                    )
           )
    fvnames(lacfv, a = ".") <- "MVL"
  }
  return(lacfv)
}


##########################
##The following function calculates lacunarity for a box with side lengths 2*bX+1 and 2*bY+1 (in pixels). The RS version is automatically calculated by ignoring those boxes that have sums that includa NA values. 
#eg mvlgb_intern.rcpproll(xiim,5,5,5*0.8)
#the obswin is only for the raw version and must be an owin object. 
#uses rcpproll
mvlgb_intern.rcpproll <- function(xiim, sidep, inclraw, obswin = Frame(xiim)){
  mat <- as.matrix(xiim)
  if ( (sidep > nrow(mat)) | (sidep > ncol(mat))){
    mvlgb.rs <- NA
    sampmean.rs <- NA
    samp2ndmom.rs <- NA
  }
  else {
    movline.overrows <- RcppRoll::roll_sum(mat, sidep)
    movline.overrowthencols <- RcppRoll::roll_sum(t(movline.overrows), sidep) * xiim$xstep * xiim$ystep
    sampmean.rs <- mean(movline.overrowthencols, na.rm = TRUE) #sample mean
    samp2ndmom.rs <- mean(movline.overrowthencols ^ 2, na.rm = TRUE) #biased sample second moment
    mvlgb.rs <- samp2ndmom.rs / (sampmean.rs ^ 2) - 1
  }

  if (inclraw){
    #lacunarity if the box centeres can be everywhere (aka no boundary correction)
    matpad <- matrix(0, ncol = ncol(mat) + 2 * sidep, nrow = nrow(mat) + 2 * sidep)
    matpad[1 * sidep + 1:nrow(xiim), 1 * sidep + 1:ncol(xiim)] <- mat
    matpad[is.na(matpad)] <- 0
    movline.overrows <- RcppRoll::roll_sum(matpad, sidep)
    movline.overrowthencols <- RcppRoll::roll_sum(t(movline.overrows), sidep)

    numpixelsinwindow <- area.owin(obswin) / (xiim$xstep * xiim$ystep)
    sampmean.raw <- sum(movline.overrowthencols) * xiim$xstep * xiim$ystep / numpixelsinwindow #note this isn't simply the mean of areafracs because the areafracs image is larger than the input image
    samp2ndmom.raw <- sum(movline.overrowthencols ^ 2) * (xiim$xstep * xiim$ystep) ^ 2 / numpixelsinwindow
    mvlgb.raw <- samp2ndmom.raw / (sampmean.raw ^ 2) - 1
  }
  if (inclraw) {return(list(raw = mvlgb.raw, RS = mvlgb.rs))}
  else {
    return(list(
      MVL = mvlgb.rs,
      s2 = samp2ndmom.rs - sampmean.rs^2,
      xbar = sampmean.rs
      ))
    }
}

#' @rdname mvlgb
lacgb <- mvlgb
