#' @title MVL estimates from shapefiles, SpatialPolygonDataFrames  and raster data
#' @description Functions for estimating MVL from regions specified in shapefiles or SpatialPolygonDataFrames given a raster file, or a \code{RasterLayer} object.
#' @export MVLests_files MVLest_region  converttologicalim  plot_MVLest_region MVLest_multipleregions 

#' @param rasterfile The path to a single layer raster file readable to RGDAL
#' @param shapefile The path to a shapefile
#' @param polydf is a SpatialPolygonsDataFrame with a single region
#' @param rasterlayer is a raster layer object
#' @param frange is list of length 2 specifying the minimum and maximum pixel values to assign as foreground. The range is inclusive.
#' @param NArange is list of length 2 specifying the minimum and maximum pixel values to assign as NA values. The range is inclusive and overrides any assignment by frange.
#' @param boxwidths The set of widths of square boxes to estimate MVL for.
#' @param estimators A list of names of MVL estimators to use. See \code{mvl()} for available list and more information
#' @param display If TRUE then \code{plot_MVLest_region} is called so that the results are automatically plotted
#' @describeIn mvlfromshapeandraster  Returns MVL estimates from regions specified in shapefile using raster data given in rasterfile
MVLests_files <- function(
  shapefile, 
  rasterfile, 
  frange, 
  NArange, 
  boxwidths, 
  estimators = c("MVLg.mattfeldt", "MVLg.pickaint",
                 "MVLcc.mattfeldt", "MVLcc.pickaint", "MVLcc.pickaH",
                 "MVLc", "MVLgb"),
  display = TRUE
){
  polysdf <- rgdal::readOGR(dirname(shapefile), sub(".shp", "", basename(shapefile)), verbose = FALSE)
  rasterlayer <- raster::raster(rasterfile)
  out <- MVLest_multipleregions(
    polysdf = polysdf,
    rasterlayer = rasterlayer, 
    frange = frange, 
    NArange = NArange, 
    boxwidths = boxwidths, 
    estimators = estimators,
    display = display
    )
  return(out)
}


#' @param polysdf A SpatialPolygonsDataFrame with each feature corresponding to a different region.
#' @describeIn mvlfromshapeandraster Splits a SpatialPolygonDataFrame into component features and estimates MVL using \code{rasterlayer} data in each of these features.
MVLest_multipleregions <- function(polysdf, 
                                       rasterlayer, 
                                       frange, 
                                       NArange, 
                                       boxwidths, 
                                       estimators = c("MVLg.mattfeldt", "MVLg.pickaint",
                                             "MVLcc.mattfeldt", "MVLcc.pickaint", "MVLcc.pickaH",
                                             "MVLc", "MVLgb"),
                                       display = TRUE){
  lpolydf <- unlistSpatialPolygonsDataframe(polysdf)
  names(lpolydf) <- as.data.frame(polysdf)[,1]
  out <- lapply(lpolydf, MVLest_region, rasterlayer = rasterlayer,
         frange = frange, NArange = NArange, boxwidths = boxwidths, estimators = estimators)
  if (display) {
    plot_MVLest_allregions(out, estname = "mvlcc.pickaH", main = "Class images and PickaH MVL estimate")
  }
  return(out)
}

#' @describeIn mvlfromshapeandraster Computes MVL estimates assuming that SpatialPolygonsDataFrame (polydf) representing a single region. 
#' Returns a list of the three items (1) the logical image used, (2) the MVL estimates and (3) the polydf used an input.
MVLest_region <- function(polydf, rasterlayer,
                          frange, NArange, boxwidths, estimators = c("MVLg.mattfeldt", "MVLg.pickaint",
                                                                      "MVLcc.mattfeldt", "MVLcc.pickaint", "MVLcc.pickaH",
                                                                      "MVLc", "MVLgb")){
  xiim <- converttologicalim(polydf, rasterlayer, frange, NArange)
  mvl.ests <- mvl(xiim, boxwidths, estimators = estimators)
  return(list(classimage = xiim, polydata = polydf, mvl.est = mvl.ests))
}

#' @describeIn mvlfromshapeandraster Converts data contained in \code{rasterlayer} into a \code{spatstat im} object with logical values of TRUE (for foreground), FALSE (for background) and NA.
converttologicalim <- function(polydf, rasterlayer, frange, NArange){
  obswin <- as.owin(polydf)
  rast <- raster::crop(rasterlayer, raster::extent(polydf)) #crop the rasterlayer to the observation window
  imobj <- as.im(rast)
  imobj[setminus.owin(Frame(imobj), obswin)] <- NA
  fground <- solutionset( (frange[[1]] <= imobj) & (imobj <= frange[[2]]) )
  NAground <- solutionset( (NArange[[1]] <= imobj) & (imobj <= NArange[[2]]) )
  
  imobj[!is.na(imobj$v)] <- FALSE #anything that isn't NA set to 0
  imobj[fground] <- TRUE
  imobj[NAground] <- NA
  return(imobj)
}

#' @describeIn mvlfromshapeandraster A special plotting function for the results of \code{MVLest_region}. Plots the image used, the MVL estimates and includes polygon attributed below plots.
#' @param returnedlist The results of a call to \code{MVLest_region}
#' @param plot.im.args A named list of arguments passed to \code{plot.im} for plotting the class image.
#' @param plot.mvl.args A names list of argument passed to \code{plot.fv} for plotting the MVL estimates.
plot_MVLest_region <- function(returnedlist, plot.im.args = list(main = "Class Image", axes = TRUE), plot.mvl.args = NULL){
  graphics::par(mfrow = c(1, 2), oma = c(2, 0, 2, 0))
  do.call(plot.im, args = c(list(x = returnedlist$classimage), plot.im.args))
  do.call(plot.fv, args = c(list(x = returnedlist$mvl.est), plot.mvl.args))
  attrchar <- paste("Region Attributes:", as.character(returnedlist$polydata@data))
  graphics::mtext(text = attrchar, side = 1, outer = TRUE, line = 0)
}

plot_MVLest_allregions <- function(returnedlist, estname = "mvlcc.pickaH", ...){
  fvall <- createfvofallregions(returnedlist, estname = estname)
  ims <- lapply(returnedlist, function(x) x$classimage)
  alist <- do.call(anylist, c(ims, list(MVL = fvall)))
  plot.anylist(alist,
               ...
               )
}

createfvofallregions <- function(output, estname = "mvlcc.pickaH"){
  fvs <- lapply(output, function(x) x$mvl.est)
  lacfv <- collapse.fv(fvs,  different = c(estname))
  return(lacfv)
}

