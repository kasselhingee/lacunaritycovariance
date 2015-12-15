library(maptools) #with rgeos installed (for ubuntu need to install GEOS separately - waiting to see if needed!)
library(raster)
library(spatstat)

#read boundary in and convert to spatstat window
boundaryMAPTOOLS <- readShapeSpatial("data/poly01_polygons.shp")
boundaryMAPTOOLS <- readShapeSpatial("../../Data/CBC_Feed_areas_requiring_investigation_split.shp")
boundaryOWIN <- as.owin(boundaryMAPTOOLS[16408,1]) #this is actually row 16407 in the QGIS viewer
boundaryOWIN <- as.owin(boundaryMAPTOOLS)  #

##read in class map convert to spatstat window
XiRASTER <- raster("data/inputBinaryMap.ers")
XiRASTER <- raster("C:/CCI-02_Work/processing_102/UM2009/2035SW_moorSE/2009_mar_UM_ucd03_2035SW_moorSE_gda94_mga50_tre_raw_CSIRO_2012092313_ver01.ers")
XiRASTER <- crop(XiRASTER,boundaryMAPTOOLS[16408,1]) #this process takes some time, but doesn't actually create anything in memory
plot(XiRASTER)
plot(boundaryOWIN,add=TRUE)
XiIMAGE <- as.im(XiRASTER) #the data is now saved into memory. It is pretty inefficient - 149MB for a 39MB ERMapper file
XiOWIN <- as.owin(XiIMAGE)

#check that area of a mask owin is number of pixels in class
area.owin(XiOWIN)
sum(XiOWIN$m)*XiOWIN$xstep*XiOWIN$ystep  #also works to use as.matrix(XiOWIN)


#estimate exact variance
varPexactEst(covprobEstimate,covarianceMap,boundaryOWIN)


maxXshiftdistance <- 3.1
c(-rev(seq(0,maxXshiftdistance,by=0.2)),seq(0.2,maxXshiftdistance,by=0.2))



plot(shift.owin(Xiinside,vec=shiftVector))
plot(boundaryOWIN,add=TRUE)
plot(rData)
plot(polyWindowowin,add=TRUE)


###############heather data
data(heather)
XiOWIN <- heather$coarse
windowOWIN <- owin(xrange=heather$coarse$xrange,yrange=heather$coarse$yrange)

coverageProb <- covpest(XiOWIN,windowOWIN)

covariancePt <- covarianceEstAtPoint(XiOWIN,windowOWIN,c(3,5))

covarianceDirectEst <- covarianceMapEst_direct(XiOWIN,windowOWIN,1,1)
filled.contour(covarianceDirectEst$covariance)


covarianceFcn <- covarianceRACS(XiOWIN,windowOWIN)
plot(covarianceFcn$covariance)
plot(covarianceFcn$covariance - coverageProb*coverageProb)

Z <- HestInPolygonalWindow(XiOWIN,windowOWIN)
plot(Z)
plot(add=TRUE,Hest(XiOWIN))#since window is rectangular should be the same


###########################
#simulating a boolean model of discs of random radius with lognormal distribution
w <- owin(xrange=c(0,10),yrange=c(0,10))
#have to simulate in a much larger area than the observation window (because grains with centres outside the window should still be observed)
lambda <- 1
pp <- rpoispp(lambda,win=w,nsim=1,drop=TRUE)

radius <- rlnorm(pp$n,meanlog=-1,sdlog=0.5) #prepare a random radius for each point
pointlocations <- cbind(X=pp$x,Y=pp$y)
pointlocations <- split(cbind(pointlocations),row(pointlocations)) #split matrix into a list of the rows
grains <- mapply(disc,radius = radius,centre=pointlocations,SIMPLIFY=FALSE) #calculate grains with their locations

#take union of all grains
xi <- union.owin(as.solist(grains))

