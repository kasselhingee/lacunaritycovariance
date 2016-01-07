library(maptools) #with rgeos installed (for ubuntu need to install GEOS separately - waiting to see if needed!)
library(raster)
library(spatstat)
#also require rgdal to be installed correctly!


#read boundary in and convert to spatstat window
boundaryMAPTOOLS <- readShapeSpatial("data/poly01_polygons.shp")
boundaryOWIN <- as.owin(boundaryMAPTOOLS)  #

##read in class map convert to spatstat window
XiRASTER <- raster("data/inputBinaryMap.ers")
XiIMAGE <- as.im(XiRASTER)
XiOWIN <- as.owin(XiIMAGE)

#check that area of a mask owin is number of pixels in class
area.owin(XiOWIN)
sum(XiOWIN$m)*XiOWIN$xstep*XiOWIN$ystep  #also works to use as.matrix(XiOWIN)


#intersect boundary with Xi observation
Xiinside = intersect.owin(XiOWIN,boundaryOWIN)   #tried using raster package directly and didn't work. Couldn't see anything in maptools
area.owin(Xiinside)
covprobEstimate = area.owin(Xiinside)/area.owin(boundaryOWIN)


#estimate covariance of a particular discrete vector v.
shiftVector = c(2,0.2) #in units of Xiinside
numerator = area.owin(intersect.owin(Xiinside,shift.owin(Xiinside,vec=shiftVector)))
denominator = area.owin(intersect.owin(boundaryOWIN,shift.owin(boundaryOWIN,vec=shiftVector)))
numerator/denominator


#experiment with plotting a map of covariance.
covarianceMapEst <- function(vecx,vecy,Xi,w){
  covAtv=covarianceEst(Xi,w,c(vecx[1],vecy[1]))
  print(c(covAtv,vecx=vecx,vecy=vecy))
  return(covAtv[3])
}
covarianceEst(XiOWIN,boundaryOWIN,shiftVector)[3]
covarianceMapEst(XiOWIN,boundaryOWIN,shiftVector[1],shiftVector[2])
shiftVectorX = 1:3*0.2
shiftVectorY = 1:3*0.2

covarMap = outer(shiftVectorX,shiftVectorY,FUN = "covarianceMapEst",XiOWIN,boundaryOWIN)


#using outer hasn't seemed to work - it needs shift.owin to accept vectorised inputs.
shiftVectorX = -30:30*0.2
shiftVectorY = -30:30*0.2
covarMap = matrix(nrow=length(shiftVectorY),ncol=length(shiftVectorX))
numeratorMap = matrix(nrow=length(shiftVectorY),ncol=length(shiftVectorX))
denominatorMap = matrix(nrow=length(shiftVectorY),ncol=length(shiftVectorX))
for (i in 1:length(shiftVectorX)){
  for (j in 1:length(shiftVectorY)){
    covarianceEstimation <- covarianceEst(XiOWIN,boundaryOWIN,c(shiftVectorX[i],shiftVectorY[j]))
    covarMap[i,j] <- covarianceEstimation["covarianceEst"]
    numeratorMap[i,j] <- covarianceEstimation["numerator"]
    denominatorMap[i,j] <- covarianceEstimation["denominator"]
  }
}
filled.contour(shiftVectorX,shiftVectorY,covarMap)

covarIm <- im(covarMap,xcol=shiftVectorX,yrow=shiftVectorY)
#trying to get a matrix of lists.
covarianceMap = list(covariance = covarMap, numerator=numeratorMap, denominator = denominatorMap, xcol = shiftVectorX, ycol = shiftVectorY, xstep=0.2,ystep=0.2)
#try array(list(NULL), c(3,2))
#or multidim array array(data = NA, dim = length(data), dimnames = NULL)
covarianceMap2 <- covarianceMapEst(XiOWIN,boundaryOWIN,6.01,6.01)
filled.contour(covarianceMap2$covariance) 

plot(shift.owin(Xiinside,vec=shiftVector))
plot(boundaryOWIN,add=TRUE)
plot(rData)
plot(polyWindowowin,add=TRUE)

#try out setcov: estimator is simply setcov(Xi)/setcov(W) #seem to be very similar
load("covarianceMapSmallBinaryMapBench.RData")
layout(matrix(c(1,2,3,4),nrow=2,ncol=2,byrow=TRUE))
numerator <- setcov(Xiinside)
plot(numerator)
image(covarianceMapSmallBinaryMapBench$numerator)
denominator <- setcov(boundaryOWIN)
plot(denominator)
image(covarianceMapSmallBinaryMapBench$denominator)
#filled.contour(covarianceMapSmallBinaryMapBench$numerator)
#filled.contour(covarianceMapSmallBinaryMapBench$denominator)
#trouble is setcov comes back with different resolutions to the setcov of the boundary!? use eval.im and "harmonize=TRUE"! to resample so they are compatible
covariance <- eval.im(numerator %/% denominator,harmonize=TRUE)#seems error prone - denominator has things like 2^-20 in it!!
harmonised <- harmonise(numerator,denominator)
plot(harmonised[[1]])
plot(harmonised[[2]])
min(denominator)
windowarea <- area.owin(boundaryOWIN)
denominatorThresh <- denominator
denominatorThresh[denominator<0.1*windowarea] <- NA
plot(denominatorThresh)

plot(eval.im(numerator / denominatorThresh,harmonize=TRUE))
image(covarianceMapSmallBinaryMapBench$covariance)

#try out function
covarianceInfo <- covarianceRACS(XiOWIN,boundaryOWIN,setCovBoundaryThresh<-0.01*area.owin(boundaryOWIN))
plot(covarianceInfo$covariance)
#compare to old funct
covarianceMap2 <- covarianceMapEst_direct(XiOWIN,boundaryOWIN,6.01,6.01)
image(covarianceMap2$covariance)

#compare to limit (assuming mixing)
asympLim <- covprobEstimate*covprobEstimate
plot(eval.im(X - asympLim,list(X=covarianceInfo$covariance)))
plot(Xiinside)
plot(add=TRUE,boundaryOWIN)


variance <- exactvarP(XiOWIN,boundaryOWIN) #came out negative for small Poly01! Probably fair!



