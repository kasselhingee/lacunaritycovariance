#edge density scratch
library(maptools) #with rgeos installed (for ubuntu need to install GEOS separately - waiting to see if needed!)
library(raster)
library(spatstat)
#also require rgdal to be installed correctly!


#read boundary in and convert to spatstat window
boundaryMAPTOOLS <- readShapeSpatial("data/poly01_polygons.shp")
boundaryOWIN <- as.owin(boundaryMAPTOOLS)  #

##read in class map convert to spatstat window
XiRASTER <- raster("data/inputBinaryMap.ers")
XiRASTER <- raster("/home/tearcor/PhDLargeDataFies/Balcatta_Park/BalcattaPark_veg_raw.ers")
XiRASTER <- raster("/home/kassel/LargeData/Balcatta_Park/BalcattaPark_veg_raw.ers")
XiIMAGE <- as.im(XiRASTER)
XiOWIN <- as.owin(XiIMAGE)

Xiinside = intersect.owin(XiOWIN,boundaryOWIN)
plot(Xiinside)

#plan - use whatever convexify uses to approximate pixels with polygons. Measure lengths of polygons.
#simulate PP on polygon boundaries too

eps=NA
#if (!is.polygonal(Xiinside)) {
  if (is.na(eps)) #if (missing(eps)) 
    eps <- diameter(Frame(Xiinside))/20
    eps <- diameter(Frame(Xiinside))/40 #this one seemed to create slightly better looking patches for inputBinaryMap
  ptm <- proc.time()
  XiinsideSimplified <- simplify.owin(Xiinside, eps) #I think only works on binary masks in dev version of spatstat #simplify.owin for balcattaPark took 3.6!! hours
  proc.time() - ptm
  balcattaParkInsideTreeSimplified <- XiinsideSimplified
  save(balcattaParkInsideTreeSimplified,file="balcattaParkInsideTreeSimplified.RData")
  #}
layout(matrix(c(1,2),nrow=1,ncol=2))
plot(Xiinside)
plot(XiinsideSimplified,color="green")
#notice the loss of a few pieces - are they important?!
warning("biased perimenter estimation!")
perimeter(XiinsideSimplified)

perimeter(as.polygonal(Xiinside))

#start extracting the edges for simulation of point process
edgelist <- edges(XiinsideSimplified,window <- boundaryOWIN) #so that edges list has the correct observation boundary (otherwise just uses smallest rectangle)
#should check that this edge list doesn't cross with itself? Don't think is necessary but could use selfcrossing.psp and/or selfcut.psp
lambda = 3
coxRealisations <- rpoisppOnLines(lambda,edgelist,nsim=6)
plot(coxRealisations)
plot(add=TRUE,XiinsideSimplified)

reduced2ndMomentMeasure <- Kmeasure(coxRealisations[[1]],diameter(Frame(boundaryOWIN))/10)
plot(reduced2ndMomentMeasure)
