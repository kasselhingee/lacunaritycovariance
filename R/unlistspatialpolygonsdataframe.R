#' @title unlist SpatialPolygonsDataFrame
#' @export unlistSpatialPolygonsDataframe
#' @description  Takes a SpatialPolygonsDataFrame object and splits it into individual SpatialPolygonsDataFrames. One for each polygon.
#' There is a nicer way: use split(, f = 1:SPDFlength).
#' 
#' @param spdf A SpatialPolygonsDataFrame
#' @return a list of SpatialPolgyonsDataFrames, each data frame contains exactly one polygon.
unlistSpatialPolygonsDataframe <- function(spdf){
  polylist <- list()
  for (i in 1:nrow(spdf)){
    polylist <- c(polylist,spdf[i,])
  }
  return(polylist)
}
