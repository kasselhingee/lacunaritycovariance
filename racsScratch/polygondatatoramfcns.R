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


