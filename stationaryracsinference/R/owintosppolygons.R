#' @title convert owin polygons to an sp SpatialPolygonsDataFrame
#' @export as.Polygon.owin  as.Polygons.Polygon   as.spdf.owin
#' 
#' @description functions for converting an owin polygon or rectangle to the sp polygon formats. This is useful for interfacing with rgdal and raster
#' 
#' @details \code{as.Polygon.owin} for converting owin objects to sp polygons without an projection information
#' 
#' @param owinpoly A single spatstat owin object of type polygonal or rectangular

#' @examples 
#' data(polygontest)
#' tess <- quadrats(maptools::as.owin.SpatialPolygons(polygontest))
#' tilelist <- tiles(tess)
#' 
#' @section Note: functions such \code{setAs()} seem much better than these current things. Also relies of owin polygons being single, closed polygons.
#' maptools has a function "as.SpatialPolygons.tess" in its source code (file spatstat1.R), but it doesn't seem to be exported.


#' @return \code{as.Polygon.owin} returns a single sp polygon object
#' 
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
  sppolygon <- rgdal::Polygon(coords,hole=FALSE)
  return(sppolygon)
}

#' @describeIn as.Polygon.owin From a single sp \code{Polygon} object returns an sp \code{Polygons} object
#' @param sppolygon is a "Polygon" object for sp.
#' @param id ids to be assigned to the polygon(s). Default is NA
#' @return an sp "Polygons" object from a single Polygon object and an id
as.Polygons.Polygon <- function(sppolygon, id="NA"){
  sppolygons <- Polygons(list(sppolygon),id)
  return(sppolygons)
}

#' @describeIn as.Polygon.owin Returns a SpatialPolgyonsDataFrame given a list of owin polygons, ids and projection information
#' @param owinpoly is a list of owin polygons
#' @param proj4string Projection information see sp and rgdal for more detail
#' @param auxdata dataframe of auxilary information for each polygon, must include an id column. Must have rows corresponding to each polygon.
#' @param idcolumn a character refering to the column to use for ids in auxdata. If none is supplied then uses rownames
#' @param match.ID see help for \code{SpatialPolygonsDataFrame}. If TRUE check that row names of auxdata match to id possibly with rearrangement
#' If character match polygon ids with a particular column in auxdata
#' 
#' @seealso SpatialPolygonsDataFrame
as.spdf.owin <- function(owinpolys, auxdata, idcolumn = NULL, proj4string=CRS(as.character(NA)), match.ID = TRUE){
  spPolys <- lapply(owinpolys,as.Polygon.owin) #convert into list of "Polygon" format
  if (is.null(idcolumn)) {id <- row.names(auxdata)}
  else {id <- auxdata$id}
  spPolys <- mapply(as.Polygons.Polygon, spPolys, id) #convert into "Polygons" format
  spPolys <- SpatialPolygons(spPolys, proj4string = proj4string)
  spdf <- SpatialPolygonsDataFrame(spPolys, auxdata, match.ID = match.ID)
  return(spdf)
}


