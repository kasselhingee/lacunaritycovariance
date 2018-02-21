#' @title Variance Estimates for Observed Area
#' @export allsae
#' @description Estimates the variance of the area of a cover type observed in a thematic map created using a fallible classifier from remote sensing.
#' @author{Kassel Hingee}

#' @param xi An observation of the RACS of interest in 1, 0, or NA valued pixels.
#' Pixels must be square. Must be a spatstat im object.
#' @param obswin Observation window
#' @param corrrad Radius of the step function in the correlation (in the same units as xi)
#' @param corrstepheight Height of the step in the correlation
#' @param n11 The number of samples that were truly class 1 and fallibly classified as class 1
#' @param n21 The number of samples that were truly class 2 and fallibly classified as class 1
#' @param n12 The number of samples that were truly class 1 and fallibly classified as class 2
#' @param n22 The number of samples that were truly class 2 and fallibly classified as class 2

#' @examples
#' xi <- heather$coarse
#' obswin <- Frame(heather$coarse)
#' corrrad <- 2
#' corrstepheight <- 0.95
#' n11 <- 90
#' n21 <- 10
#' n22 <- 95
#' n12 <- 5
#' estimates <- allsae(xi, obswin, corrrad, corrstepheight, n11, n21, n12, n22)
#' \dontrun{
#' plot.new()
#' plot.window(ylim = c(min(estimates$areahat - 2*sqrt(estimates$varhat)),
#'                      max(estimates$areahat + 2*sqrt(estimates$varhat))),
#'            xlim = c(0, length(estimates$areahat) + 1))
#' axis(1, at = 1:length(estimates$areahat), labels = row.names(estimates))
#' axis(2)
#' arrows(1:length(estimates$areahat),
#'   y0 = estimates$areahat - 2*sqrt(estimates$varhat),
#'   y1 = estimates$areahat + 2*sqrt(estimates$varhat),
#'   length = 0.05, angle = 90, code =3)
#' }


allsae <- function(xi, obswin, corrrad, corrstepheight, erosionrad, n11, n21, n12, n22){
  #probability that a labelled tree is not really tree
  p21 <- n21 / (n11 + n21)
  #probability that a labelled non-tree is reall tree
  p12 <- n12 / (n22 + n12)

  results <- data.frame(NULL)
  results["Confusion Matrix Only", "areahat"] <- 
    sae.v1ab.mean(area.owin(xi), area.owin(obswin), n11, n21, n12, n22)
  results["Confusion Matrix Only", "varhat"] <-
    sae.v1ab.var(area.owin(xi), area.owin(obswin), n11, n21, n12, n22)

  results["IID Pixels", "areahat"] <- 
    sae.v1b.mean(area.owin(xi), area.owin(obswin), xi$xstep, p21, p12)
  results["IID Pixels", "varhat"] <- 
    sae.v1b.var(area.owin(xi), area.owin(obswin), xi$xstep, p21, p12)
  
  results["IID Pixels + C Mat", "areahat"] <- 
    sae.v1b.wsu.mean(area.owin(xi), area.owin(obswin), xi$xstep, n11, n21, n12, n22)
  results["IID Pixels + C Mat", "varhat"] <- 
    sae.v1b.wsu.var(area.owin(xi), area.owin(obswin), xi$xstep, n11, n21, n12, n22)

  results["IID Regions", "areahat"] <-
    sae.v2d.mean(area.owin(xi), area.owin(obswin), xi$xstep, p21, p12, corrrad)
  results["IID Regions", "varhat"] <-
    sae.v2d.var(area.owin(xi), area.owin(obswin), xi$xstep, p21, p12, corrrad)
  
  results["IID Regions + C Mat", "areahat"] <-
    sae.v2d.wsu.mean(area.owin(xi), area.owin(obswin), xi$xstep, n11, n21, n12, n22, corrrad)
  results["IID Regions + C Mat", "varhat"] <-
    sae.v2d.wsu.var(area.owin(xi), area.owin(obswin), xi$xstep, n11, n21, n12, n22, corrrad)

  results["Step Covariance","areahat"] <- sae.v3.mean(xi, obswin, p21 = p21, p12 = p12)
  results["Step Covariance", "varhat"] <- sae.v3.var(xi, obswin, corrrad, corrstepheight, p21 = p21, p12 = p12)

  results["Step Cov + C Mat","areahat"] <- sae.v3.wsu.mean(xi, obswin, n11, n21, n12, n22)
  results["Step Cov + C Mat", "varhat"] <- sae.v3.wsu.var(xi, obswin, corrrad, corrstepheight, n11, n21, n12, n22)

  results["ErodeDilate", "areahat"] <- sae.v4.mean(xi, obswin)
  results["ErodeDilate", "varhat"] <- sae.v4.var(xi, obswin, erosionrad)
  
  return(results)
}
