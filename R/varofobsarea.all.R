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


allsae <- function(xi, obswin, corrrad, corrstepheight, n11, n21, n12, n22){
  #probability that a labelled tree is not really tree
  p21 <- n21 / (n11 + n21)
  #probability that a labelled non-tree is reall tree
  p12 <- n12 / (n22 + n12)

  results <- data.frame(NULL)
  results["sae.v1ab", "areahat"] <- 
    sae.v1ab.mean(area.owin(xi), area.owin(obswin), n11, n21, n12, n22)
  results["sae.v1ab", "varhat"] <-
    sae.v1ab.var(area.owin(xi), area.owin(obswin), n11, n21, n12, n22)

  results["sae.v1b", "areahat"] <- 
    sae.v1b.mean(area.owin(xi), area.owin(obswin), xi$xstep, p21, p12)
  results["sae.v1b", "varhat"] <- 
    sae.v1b.var(area.owin(xi), area.owin(obswin), xi$xstep, p21, p12)
  
  results["sae.v1bb", "areahat"] <- 
    sae.v1b.wsu.mean(area.owin(xi), area.owin(obswin), xi$xstep, n11, n21, n12, n22)
  results["sae.v1bb", "varhat"] <- 
    sae.v1b.wsu.var(area.owin(xi), area.owin(obswin), xi$xstep, n11, n21, n12, n22)

  results["sae.v2d", "areahat"] <-
    sae.v2d.mean(area.owin(xi), area.owin(obswin), xi$xstep, p21, p12, corrrad)
  results["sae.v2d", "varhat"] <-
    sae.v2d.var(area.owin(xi), area.owin(obswin), xi$xstep, p21, p12, corrrad)

  results["sae.v3","areahat"] <- sae.v3.mean(
             as.im(xi, na.replace = 0)[obswin, drop = FALSE],
             obswin,
             p21 = p21,
             p12 = p12)
  results["sae.v3", "varhat"] <- sae.v3.var(
                as.im(xi, na.replace = 0)[obswin, drop = FALSE],
                obswin = obswin,
                corrrad = corrrad,
                corrstepheight = corrstepheight,
                p21 = p21,
                p12 = p12)
  
  return(results)
}
