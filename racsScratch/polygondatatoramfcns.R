#accept list of polygons, all in one file. *Efficiently* gets their values in memory eventually

#' @depends raster
#' @export putencompassingrastervaluesinram

#' @param polygons is a SpatailPolygonsDataFrame
#' @param inrasterfilenam The file to extract raster values from. May only work for single band files
getrastervaluesofpolygons <- function(polygons,inrasterfilename){
  encompassingraster <- putencompassingrastervaluesinram( polygons,inrasterfilename)
  
  polylist <- unlistSpatialPolygonsDataframe(polygons)
  #extract raster values for each polygon:
  polygonrasters= mapply(crop,polylist,MoreArgs=list(x=encompassingraster),simplify=FALSE)
  
  return(polygonrasters)
}




#' @param polygons is a SpatialPolygonDataFrame.
#' @param filename is the location of a raster file that contains *all* of the polygons
putencompassingrastervaluesinram <- function(polygons,filename){
  stopifnot(class(polygons)[1]=="SpatialPolygonsDataFrame")
  rasterobject <- raster(inrasterfilename)
  if (!(extent(rasterobject) >= extent(polygons))){
    stop("Raster file doesn't cover all polygons. Exiting")
  }

  rasterinram <- transfertoram(rasterobject,extent(polygons),2)
  return(rasterinram)
}

#' @param raster Any raster object.
#' @param sizefactor A multiplicative safety margin - it must be possible to hold \code{sizefactor}*raster data size in RAM.
#' @param extent Giving the extent of values to load into ram
transfertoram <- function(raster,extent, sizefactor){
  rasterforextraction <- crop(raster,extent)
  #crop is really slow!!
  if (sizefactor < 2){sizefactor <- 2} #for the current version of this function need to be able to put two copies into RAM
  small <- canProcessInMemory(rasterforextraction,2)
  if (!small) {stop("raster too large to copy to RAM. Exiting")}
  values <- getValues(rasterforextraction,format="matrix")
  rasterinram <- raster(values,template=rasterforextraction)
  rm(values) #two copies in RAM now!
  return(rasterinram)
  
  #It should be possible to improve this function by using Rgdal directly.
}

#' @details Takes a SpatialPolygonsDataFrame object and splits it into individual SpatialPolygonsDataFrames. One for each polygon.
#' I wish there was a nicer way than currently in this function.
unlistSpatialPolygonsDataframe <- function(spdf){
  polylist <- list()
  for (i in 1:nrow(spdf)){
    polylist <- c(polylist,spdf[i,])
  }
  return(polylist)
}


######################################
#' @description Functions for taking a list of owin polygons/rectangles and converting to a SpatialPolygonsDataFrame

#' @param x A spatstat owin object of type polygonal or rectangular
#only works for polygonal and rectangular owin
#this particular heirarchy of classes from sp make no sense to me!
as.Polygon.owin <- function(x){
  if (is.polygonal(x)){
    coords <- matrix(c(x$bdry[[1]]$x,x$bdry[[1]]$y),byrow=FALSE,nrow=length(x$bdry[[1]]$x),ncol=2)
    coords <- rbind(coords,coords[1,]) #finish loop as per requirements of Polygon function
  }
  else if (is.rectangle(x)){
    coords <- matrix(c(
      x$xrange[1],x$yrange[1],
      x$xrange[1],x$yrange[2],
      x$xrange[2],x$yrange[2],
      x$xrange[2],x$yrange[1],
      x$xrange[1],x$yrange[1]
    ),byrow=TRUE,ncol=2,nrow=5)
  }
  else {
    stop(paste("Error: can only convert polygon or rectangle type owin's to spatial Polygons."))
  }
  sppolygon <- Polygon(coords,hole=FALSE)
  return(sppolygon)
}

#' @param x is a "Polygon" object for sp.
#' @param id an id to be assigned to the polygon. Default is NA
#' @param returns a "Polygons" object
as.Polygons.Polygon <- function(x, id="NA"){
  sppolygons <- Polygons(list(x),id)
}
