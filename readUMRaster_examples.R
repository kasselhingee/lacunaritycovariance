#given a polygon, function works out which UM tile is needed and returns raster data

require(raster)
require(maptools)
require(rgdal)
#scratch stuff
#polygon with same projection as UM data
polyOGR <- readOGR("data","poly03") #works if GDA94 / MGA 50 are selected as the CRS when creating layer in QGIS
crs(proj4string(polyOGR))
polygonMAPTOOLS <- readShapeSpatial("data/poly03.shp",proj4string = crs(proj4string(polyOGR)))

#same polygon with different projection to UM data
polyOGR <- readOGR("data","poly03_albers") #in CRS GDA94 / australian albers
crs(proj4string(polyOGR))
polygonMAPTOOLS <- readShapeSpatial("data/poly03_albers.shp",proj4string = crs(proj4string(polyOGR)))


proj4string(polygonMAPTOOLS) == "+proj=utm +zone=50 +south +ellps=GRS80 +units=m +no_defs"


mspImage <- readUMraster(polygonMAPTOOLS,"dom","C:/CCI-02_Work/processing_102/UM2009/")
treMask <- readUMraster(polygonMAPTOOLS,"tre","C:/CCI-02_Work/processing_102/UM2009/")

plot(mspImage)

"+proj=utm +zone=50 +south +ellps=GRS80 +units=m +no_defs"

plot(add=TRUE,spTransform(polygonMAPTOOLS,CRS("+proj=utm +zone=50 +south +ellps=GRS80 +units=m +no_defs")))

#####Numerical test of spTransform
layout(matrix(c(1,2),ncol=2))




polygon_albers == polygon_MGA50
plot(polyOGR_MGA50,col ="red")
plot(add=TRUE,polygon_albers_trans,col ="green")

proj4string(polygonMAPTOOLS) == "+proj=utm +zone=50 +south +ellps=GRS80 +units=m +no_defs"

#test reprojection:
test_reprojection()

