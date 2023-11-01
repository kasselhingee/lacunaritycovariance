#' @keywords internal
#' @details 
#' Random closed sets (RACS) (Chiu et al., 2013; Molchanov, 2005) are a well known tool for modelling binary coverage maps. 
#' The package author recently developed new, improved estimators of gliding box lacunarity (GBL) for RACS (Hingee et al., 2017) and described contagion-like properties for RACS (Hingee, 2016).
#' The PhD thesis (Hingee, 2019) provides additional background for GBL, and for RACS in landscape metrics (which includes contagion).
#'
#' This package expects RACS observations to be in the form of binary maps either in raster format, or as a set representing foreground with a second set giving the observation window.
#' If in raster format, the binary map is expected to be a [spatstat.geom::im()] object with pixel values that are only 1 and 0, or are logically valued (i.e. `TRUE` or `FALSE`). In both cases the observation window is taken to be the set of pixels with values that are not `NA` (i.e. `NA` values are considered outside the observation window).
#' The foreground of the binary map, corresponding to locations within the realisation of the RACS, is taken to be pixels that have value `1` or `TRUE`.
#' If the binary map is in set format then a [spatstat.geom::owin()] object is used to represent foreground and a second `owin` object is used to represent the observation window.
#'
#' We will usually denote a RACS as \eqn{\Xi} ('Xi') and a realisation of \eqn{\Xi} observed as a binary map as \eqn{xi}. We will usually denote the observation window as `obswin`.
#'
#' A demonstration converting remotely sensed data into a binary map in `im` format can be accessed by typing `demo("import_remote_sense_data", package = "lacunaritycovariance")`.
#' A short example of estimating RACS properties can be found in the vignette `estimate_RACS_properties`, which can be accessed with `vignette("estimate_RACS_properties")`.
#'
#' The key functions within this package for estimating properties of RACS are:
#' * [coverageprob()] estimates the coverage probability of a stationary RACS
#' * [racscovariance()] estimates the covariance of a stationary RACS
#' * [gbl()] estimates the GBL of a stationary RACS
#' * [cencovariance()] estimates the centred covariance of a stationary RACS
#' * [paircorr()] estimates the pair-correlation of a stationary RACS
#' * [secondorderprops()] estimates GBL, covariance and other second order properties of stationary RACS
#' * [contagdiscstate()] estimates the disc-state contagion of a stationary RACS
#'
#' Key functions for simulating RACS are:
#' * [rbdd()] simulates a Boolean model with grains that are discs with fixed radius (deterministic discs).
#' * [rbdr()] simulates a Boolean model with grains that are rectangles of fixed size and orientation.
#' * [rbpto()] simulates a Boolean model with grains that of fixed shape and random scale distributed according to a truncated Pareto distribution.
#' * [placegrainsfromlib()] randomly places grains on a set of points (used to simulate Boolean models and other germ-grain models).
#' @references
#' Chiu, S.N., Stoyan, D., Kendall, W.S. and Mecke, J. (2013) *Stochastic Geometry and Its Applications*, 3rd ed. Chichester, United Kingdom: John Wiley & Sons.
#' 
#' Hingee, K.L. (2016) Statistics for patch observations. *International Archives of the Photogrammetry, Remote Sensing and Spatial Information Sciences* pp. 235-242. Prague: ISPRS.
#'
#' Hingee, K.L. (2019) *Spatial Statistics of Random Closed Sets for Earth Observations*. PhD: Perth, Western Australia: University of Western Australia. \doi{10.26182/5dbb81b6480f9}
#' 
#' Hingee K, Baddeley A, Caccetta P, Nair G (2019). Computation of lacunarity from covariance of spatial binary maps. *Journal of Agricultural, Biological and Environmental Statistics*, 24, 264-288. \doi{10.1007/s13253-019-00351-9}.
#' 
#' Molchanov, I.S. (2005) *Theory of Random Sets*. USA: Springer.
#' @keywords package
#' @keywords spatial
#' @examples
#' # Estimates from the heather data in spatstat
#' xi_owin <- heather$coarse
#' xi_owin_obswin <- Frame(heather$coarse)
#' 
#' # Convert binary map to an im object (optional)
#' xi <- as.im(xi_owin, value = TRUE, na.replace = FALSE)
#' 
#' # Estimate coverage probability, covariance, GBL, and disc-state contagion
#' cphat <- coverageprob(xi)
#' cvchat <- racscovariance(xi, estimator = "pickaH")
#' \donttest{
#'   gblhat <- gbl(xi, seq(0.1, 5, by = 1), estimators = "GBLcc.pickaH")
#'   contagds <- contagdiscstate(Hest(xi), Hest(!xi), p = cphat)
#' }
#' 
#' # Simulate a Boolean model with grains that are discs of fixed radius:
#' \donttest{
#'   xi_sim <- rbdd(10, 0.1, owin())
#' }
"_PACKAGE"


