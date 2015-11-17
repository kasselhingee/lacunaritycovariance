#test functions according to previously found values of various quanitities - great for updating code
#testing with tiny real BinaryMap and a weird shaped polygon
library(maptools) #with rgeos installed (for ubuntu need to install GEOS separately - waiting to see if needed!)
library(raster)
library(spatstat)

#read boundary in and convert to spatstat window
boundaryMAPTOOLS <- readShapeSpatial("data/poly01_polygons.shp")
boundaryOWIN <- as.owin(boundaryMAPTOOLS)  #

##read in class map convert to spatstat window
XiRASTER <- raster("data/inputBinaryMap.ers")
XiIMAGE <- as.im(XiRASTER)
XiOWIN <- as.owin(XiIMAGE)

testResultsSummary <- data.frame(success = NA, comments = "", row.names = c("coverage prob"), check.names=TRUE, stringsAsFactors = FALSE)

#test coverage probility estimate
source("covpest.R")
pest <- covpest(XiOWIN,boundaryOWIN)
testResultsSummary["coverage prob","success"] <- abs(pest - 0.332249)< 5e-8
testResultsSummary["coverage prob","comments"] <- paste("coverage estimate was ", abs(pest - 0.332249), " from benchmark value.")


#test Hest in polygonal window
source("HestInPolygonalWindow.R")
load("ZsmallBinaryMapBenchmark.RData")
Z <- HestInPolygonalWindow(XiOWIN,boundaryOWIN)
Z == ZsmallBinaryMapBenchmark

testResultsSummary
