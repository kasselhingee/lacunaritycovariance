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

