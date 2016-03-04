#' @title Functions for reading loading remote sensing data into RAM for quick extraction later
#' @export getrastervaluesofpolygons putencompassingrastervaluesinram transfertoram
#' 
#' @description When extracting multiple polygons of data it is much faster to read in a large chunk of data to ram and then extract the polygons from this chunk then to extract each polygon separately from the data on HDD. 
#' These functions are designed to help with this. Currently using raster functions, I think the best thing will be to shift to rgdal read functions in the future.
#' 

#' @param polygons is a SpatialPolygonDataFrame.
#' @param inrasterfilename is the location of a raster file that contains *all* of the polygons
getrastervaluesofpolygons <- function(polygons,inrasterfilename){
  encompassingraster <- putencompassingrastervaluesinram( polygons,inrasterfilename)
  
  polylist <- unlistSpatialPolygonsDataframe(polygons)
  #extract raster values for each polygon:
  polygonrasters= mapply(raster::crop,polylist,MoreArgs=list(x=encompassingraster),simplify=FALSE)
  
  return(polygonrasters)
}

#' @describeIn getrastervaluesofpolygons 
putencompassingrastervaluesinram <- function(polygons,inrasterfilename){
  stopifnot(class(polygons)[1]=="SpatialPolygonsDataFrame")
  rasterobject <- raster::raster(inrasterfilename)
  if (!(raster::extent(rasterobject) >= raster::extent(polygons))){
    stop("Raster file doesn't cover all polygons. Exiting")
  }
  
  rasterinram <- transfertoram(rasterobject,raster::extent(polygons),2)
  return(rasterinram)
}


#' @describeIn getrastervaluesofpolygons Tests if a chunk of data can be fit in ram, and if so transfers it to ram.
#' @param raster Any raster object.
#' @param sizefactor A multiplicative safety margin - it must be possible to hold \code{sizefactor}*raster data size in RAM.
#' @param extent Giving the extent of values to load into ram
transfertoram <- function(raster, extent, sizefactor){
  rasterforextraction <- raster::crop(raster,extent)
  #crop is really slow!!
  if (sizefactor < 2){sizefactor <- 2} #for the current version of this function need to be able to put two copies into RAM
  small <- raster::canProcessInMemory(rasterforextraction,2)
  if (!small) {stop("raster too large to copy to RAM. Exiting")}
  values <- raster::getValues(rasterforextraction,format="matrix")
  rasterinram <- raster::raster(values,template=rasterforextraction)
  rm(values) #two copies in RAM now!
  gc() #forces R to deallocate removed copy
  return(rasterinram)
  
  #It should be possible to improve this function by using Rgdal directly.
}

