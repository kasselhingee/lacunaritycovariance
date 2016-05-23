#demo read UM files


library(raster)
library(spatstat)
library(rgdal)
library(readUMdata)


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
treMask <- readUMraster(polyOGR,"tre","C:/CCI-02_Work/processing_102/UM2009/")
plot(treMask)
plot(add=TRUE, polyOGR)

covpest(treMask,w)


