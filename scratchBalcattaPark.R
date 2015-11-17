library(maptools) #with rgeos installed (for ubuntu need to install GEOS separately - waiting to see if needed!)
library(raster)
library(spatstat)

#read boundary in and convert to spatstat window
boundaryMAPTOOLS <- readShapeSpatial("C:/Users/hin117/Documents/Balcatta_Park/park_polygons.shp")
boundaryOWIN <- as.owin(boundaryMAPTOOLS)  #

##read in class map convert to spatstat window
XiRASTER <- raster("C:/Users/hin117/Documents/Balcatta_Park/BalcattaPark_veg_raw.ers")
XiIMAGE <- as.im(XiRASTER)
XiOWIN <- as.owin(XiIMAGE)

#check that area of a mask owin is number of pixels in class
area.owin(XiOWIN)
sum(XiOWIN$m)*XiOWIN$xstep*XiOWIN$ystep  #also works to use as.matrix(XiOWIN)


#intersect boundary with Xi observation
Xiinside = intersect.owin(XiOWIN,boundaryOWIN)   #tried using raster package directly and didn't work. Couldn't see anything in maptools
area.owin(Xiinside)
coverageEstimate = area.owin(Xiinside)/area.owin(boundaryOWIN)


#estimate covariance of a particular discrete vector v.
shiftVector = c(5,0.2) #in units of Xiinside
numerator = area.owin(intersect.owin(Xiinside,shift.owin(Xiinside,vec=shiftVector)))
denominator = area.owin(intersect.owin(boundaryOWIN,shift.owin(boundaryOWIN,vec=shiftVector)))
numerator/denominator



#calculate covariance using a for loop vectors up +- 6m
Sys.time()
elapsedAtStart <- proc.time()
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
proc.time() - elapsedAtStart
Sys.time()
#took time: 
##  user     system elapsed 
##  2157.87  304.97 2468.30 
##  35min    5min   41.1min
filled.contour(shiftVectorX,shiftVectorY,covarMap)
filled.contour(shiftVectorX,shiftVectorY,covarMap,zlim=c(1.5*coverageEstimate^2,max(covarMap)))


covarianceMap = list(covariance = covarMap, numerator=numeratorMap, denominator = denominatorMap)
1/(area.owin(boundaryOWIN)^2)*0.2*0.2*sum(covarianceMap$numerator-covarianceMap$denominator*covpest(XiOWIN,boundaryOWIN)^2,na.rm=TRUE)

exactVariance = varPexactEst(covpest(XiOWIN,boundaryOWIN),covarianceMap,boundaryOWIN)

