#' @title convert owin polygons to an sp SpatialPolygonsDataFrame
#' @export owin2Polygons as.SpatialPolygons.tess as.SpatialPolygons.owin    as.Polygon.owin  as.Polygons.Polygon  as.spdf.owin 
#' 
#' @description functions for converting an owin polygon or rectangle to the sp polygon formats. This is useful for interfacing with rgdal and raster
#' 
#' @details \code{as.Polygon.owin} for converting owin objects to sp polygons without an projection information
#' 

#' @examples 
#' data(polygontest)
#' tess <- quadrats(maptools::as.owin.SpatialPolygons(polygontest))
#' tilelist <- tiles(tess)
#' id <- data.frame(id = 1:length(tilelist),stringsAsFactors = FALSE)
#' tileauxdata <- data.frame(t(simplify2array(strsplit(names(tilelist), "[ , ]"))),stringsAsFactors=FALSE)
#' tileauxdata <- tileauxdata[,c(3,6)]
#' names(tileauxdata) <- c("TileRow", "TileCol")
#' tileauxdata$TileRow <- factor(tileauxdata$TileRow,ordered=TRUE) 
#' tileauxdata$TileCol <- factor(tileauxdata$TileCol,ordered=TRUE)
#' auxdata <- cbind(id,tileauxdata)
#' 
#' spdf1 <- as.spdf.owin(tilelist, tileauxdata, 
#'                      proj4string = sp::CRS(sp::proj4string(polygontest)), match.ID=FALSE)
#'                      
#' spdf2 <- sp::SpatialPolygonsDataFrame(as.SpatialPolygons.tess(tess), auxdata, match.ID = FALSE)
#' sp::proj4string(spdf2) <- sp::CRS(sp::proj4string(polygontest))
#' 
#' @section Note: functions such \code{setAs()} seem much better than these current things. Also relies of owin polygons being single, closed polygons.
#' maptools has a function "as.SpatialPolygons.tess" in its source code (file spatstat1.R), but it doesn't seem to be exported (ie. I can't use it). **Try to work out why
#' @section ToDo: test if it works on owin regions with internal polygon boundaries

#' @return \code{owin2Polygons} returns a single sp Polygons object (each owin might containg multiple polygonal boundaries)
#' @param x an owin object, or a spatstat tesselation object
#' @param id id correspond to polygons.
# methods for coercion to Spatial Polygons by Adrian Baddeley

owin2Polygons <- function(x, id="1") {
  if (!requireNamespace("spatstat", quietly = TRUE))
    stop("package spatstat required for as.owin.SpatialPixelsDataFrame")
  stopifnot(spatstat::is.owin(x))
  x <- spatstat::as.polygonal(x)
  closering <- function(df) { df[c(seq(nrow(df)), 1), ] }
  if (x$type == "polygonal") {
    pieces <- lapply(x$bdry,
                     function(p) {
                       sp::Polygon(coords=closering(cbind(p$x,p$y)),
                               hole=spatstat::is.hole.xypolygon(p))  })
  } else if (x$type == "rectangle") {
    rectCrds <- cbind(x$xrange[c(1,1,2,2,1)], x$yrange[c(1,2,2,1,1)])
    pieces <- list(sp::Polygon(rectCrds, hole=FALSE))
  } else stop("owin2Polygons: unknown type:", x$type)
  z <- sp::Polygons(pieces, id)
  return(z)
}

#' @describeIn owin2Polygons converts a tesselation to a list of Polygons
as.SpatialPolygons.tess <- function(x) {
  #require(spatstat)
  if (!requireNamespace("spatstat", quietly = TRUE))
    stop("package spatstat required for as.owin.SpatialPixelsDataFrame")
  stopifnot(spatstat::is.tess(x))
  y <- spatstat::tiles(x)
  nam <- names(y)
  z <- list()
  for(i in seq(y)) {
    zi <- try(owin2Polygons(y[[i]], nam[i]), silent=TRUE)
    if (class(zi) == "try-error") {
      warning(paste("tile", i, "defective\n", as.character(zi)))
    } else {
      z[[i]] <- zi
    }
  }
  return(sp::SpatialPolygons(z))
}

#' @describeIn owin2Polygons a wrapper of as.Polygons.owin that returns a "SpatialPolygons" object, but without any projection information
as.SpatialPolygons.owin <- function(x) {
  if (!requireNamespace("spatstat", quietly = TRUE))
    stop("package spatstat required for as.owin.SpatialPixelsDataFrame")
  #require(spatstat)
  stopifnot(spatstat::is.owin(x))
  y <- owin2Polygons(x)
  z <- sp::SpatialPolygons(list(y))
  return(z)
}




#' @describeIn owin2Polygons Kassel Hingee's original attempt to convert owins to Polygons - wont work for owin regions with boundaries that are multiple polygons
#' @param owinpoly is a polygonal or rectangular owin
as.Polygon.owin <- function(owinpoly){
  if (is.polygonal(owinpoly)){
    coords <- matrix(c(owinpoly$bdry[[1]]$x,owinpoly$bdry[[1]]$y),byrow=FALSE,nrow=length(owinpoly$bdry[[1]]$x),ncol=2)
    coords <- rbind(coords,coords[1,]) #finish loop as per requirements of Polygon function
  }
  else if (is.rectangle(owinpoly)){
    coords <- matrix(c(
      owinpoly$xrange[1],owinpoly$yrange[1],
      owinpoly$xrange[1],owinpoly$yrange[2],
      owinpoly$xrange[2],owinpoly$yrange[2],
      owinpoly$xrange[2],owinpoly$yrange[1],
      owinpoly$xrange[1],owinpoly$yrange[1]
    ),byrow=TRUE,ncol=2,nrow=5)
  }
  else {
    stop(paste("Error: can only convert polygon or rectangle type owin's to spatial Polygons."))
  }
  sppolygon <- sp::Polygon(coords,hole=FALSE)
  return(sppolygon)
}

#' @describeIn owin2Polygons From a single sp \code{Polygon} object returns an sp \code{Polygons} object
#' @param sppolygon is a "Polygon" object for sp.
#' @return an sp "Polygons" object from a single Polygon object and an id
as.Polygons.Polygon <- function(sppolygon, id="NA"){
  sppolygons <- sp::Polygons(list(sppolygon),id)
  return(sppolygons)
}

#' @describeIn owin2Polygons Returns a SpatialPolgyonsDataFrame given a list of owin polygons, ids and projection information
#' @param owinpolys is a list of owin polygons
#' @param proj4string Projection information see sp and rgdal for more detail
#' @param auxdata dataframe of auxilary information for each polygon, must include an id column. Must have rows corresponding to each polygon.
#' @param idcolumn a character refering to the column to use for ids in auxdata. If none is supplied then uses rownames
#' @param match.ID see help for \code{SpatialPolygonsDataFrame}. If TRUE check that row names of auxdata match to id possibly with rearrangement
#' If character match polygon ids with a particular column in auxdata
#' 
#' @seealso SpatialPolygonsDataFrame
as.spdf.owin <- function(owinpolys, auxdata, idcolumn = NULL, proj4string=sp::CRS(as.character(NA)), match.ID = TRUE){
  spPolys <- lapply(owinpolys,as.Polygon.owin) #convert into list of "Polygon" format
  if (is.null(idcolumn)) {id <- row.names(auxdata)}
  else {id <- auxdata$id}
  spPolys <- mapply(as.Polygons.Polygon, spPolys, id) #convert into "Polygons" format
  spPolys <- sp::SpatialPolygons(spPolys, proj4string = proj4string)
  spdf <- sp::SpatialPolygonsDataFrame(spPolys, auxdata, match.ID = match.ID)
  return(spdf)
}


