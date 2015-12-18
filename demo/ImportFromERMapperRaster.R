#demo - import from ERMapper raster files and ESRI shapefiles
#need rgdal installed correctly

library(maptools) #additional abilities of maptools useful for when shapefile doesn't include projection information
library(raster)
library(spatstat)
library(rgdal)

#######################import polygons######################
#purely with rgdal
polyOGR <- readOGR("data","poly03") #for reading the projection information AND the polygons (maptools readShapeSpatial only reads the polygons)
crs(proj4string(polyOGR))
w <- as.owin(polyOGR)


#OR using both maptools and rgdal
crsstring <- OGRSpatialRef("data","poly03")
polygonMAPTOOLS <- readShapeSpatial("data/testshp.shp",proj4string = crs(crsstring))
w <- as.owin(polygonMAPTOOLS)

#OR without any projection information
polygonMAPTOOLS <- readShapeSpatial("data/poly03.shp",proj4string = crs(crsstring))
w <- as.owin(polygonMAPTOOLS)

###########################raster################################
#import Raster Data: read in class map convert to spatstat window
XiRASTER <- raster("/home/tearcor/PhDLargeDataFies/Balcatta_Park/BalcattaPark_veg_raw.ers")
XiRASTER <- crop(XiRASTER,extent(polyOGR))
XiIMAGE <- as.im(XiRASTER)
XiOWIN <- as.owin(XiIMAGE)

#estimate RACS properties using observation:
covar <- covarianceRACS(XiOWIN,w)
plot(covar)


