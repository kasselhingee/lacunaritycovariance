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
Xiinside = intersect.owin(XiOWIN,boundaryOWIN)

testResultsSummary <- data.frame(success = NA, comments = "", row.names = c("coverage prob"), check.names=TRUE, stringsAsFactors = FALSE)

#test coverage probility estimate
source("covpest.R")
pest <- covpest(XiOWIN,boundaryOWIN)
testResultsSummary["coverage prob","success"] <- abs(pest - 0.332249)< 5e-8
testResultsSummary["coverage prob","comments"] <- paste("coverage estimate was ", abs(pest - 0.332249), " from benchmark value.")

  #test covarianceEstimate
  source("covEst.R")
  #comparing against
  #shiftVector = c(2,0.2)
  #numerator = area.owin(intersect.owin(Xiinside,shift.owin(Xiinside,vec=shiftVector)))
  #denominator = area.owin(intersect.owin(boundaryOWIN,shift.owin(boundaryOWIN,vec=shiftVector)))
  #covarianceSmallBinaryMapBench <- numerator/denominator
load("covarianceSmallBinaryMapBench.RData")
covarEst <- covarianceEst(Xiinside,boundaryOWIN,c(2,0.2))
testResultsSummary["covarianceEst","success"] <- abs(covarEst$covarianceEst - covarianceSmallBinaryMapBench)< 5e-8
testResultsSummary["covarianceEst","comments"] <- paste("covariance estimate was ", abs(covarEst$covarianceEst - covarianceSmallBinaryMapBench), " from benchmark value.")

#test covarianceMapEstimate
load("covarianceMapSmallBinaryMapBench.RData")
success <- FALSE
 covarMap <- covarianceMapEst(XiOWIN,boundaryOWIN,6.001,6.001)
 max(abs(covarMap$covariance - covarianceMapSmallBinaryMapBench$covariance),na.rm=TRUE)
 sucess <- max(abs(covarMap$covariance - covarianceMapSmallBinaryMapBench$covariance),na.rm=TRUE) < 5e-8
 testResultsSummary["covarianceMapEst","success"] <- success
 testResultsSummary["covarianceMapEst","comments"] <- paste("")
 
 
#test Hest in polygonal window
source("HestInPolygonalWindow.R")
load("HestSmallBinaryMapBenchmark.RData")
Z <- HestInPolygonalWindow(XiOWIN,boundaryOWIN)
Z == HestSmallBinaryMapBenchmark

testResultsSummary

