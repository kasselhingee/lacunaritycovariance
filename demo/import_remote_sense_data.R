
opa <- par(mfrow=c(1,1))

##################
# IMPORTING REMOTE SENSING DATA DEMO
# By Kassel Liam Hingee
##################

# This document demonstrates how to convert raster data stored in a
# remote sensing format and an observation window in an ESRI shapefile
# format into a format suitable for the package lacunaritycovariance.

# For this demo you will need the following R packages installed: 
#  rgdal You will need to install GDAL (http://www.gdal.org/)
#        first which is available for Windows, Mac and Linux.
#        An easy installation of GDAL on Windows is through
#        installing OSGeo4W (http://trac.osgeo.org/osgeo4w/).
#  maptools
#  raster
#  spatstat

# The goal is to have 
# 1) an observation window in spatstat's owin format
# 2) the raster data as a spatstat im object containing values of 1,
#    0 or NA representing foreground, background and outside the
#    observation window respectively.
# The result will be something like this:
load(system.file("extdata/egbinarymap.RData", package="lacunaritycovariance")) 
plot(egbinarymap, col = c("grey", "black"), main = "The Final Result of This Demo")
rm(egbinarymap)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# 1. Reading ESRI Shapefile in SpatialPolygonsDataFrame #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
library("rgdal") # rgdal is used here to read the ESRI shapefile into a SpatialPolygonsDataFrame
regionfilepath <- system.file("extdata", package="lacunaritycovariance")
obspoly <- readOGR(regionfilepath, "aregionofinterest", verbose = FALSE)
#For ESRI Shapefiles:
#   the first argument of readOGR is the directory containing the shape files
#   the second argument of readOGR is the filename without extension
#print the coordinate projection of the polygon data for sanity
plot(obspoly, main = "Observation Window as a SpatialPolygonsDataFrame")


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# 2.Converting Observation Window as SpatialPolygonsDataFrame into owin Format  #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
library("maptools") # maptools is used to convert a SpatialPolygonsDataFrame into an owin object
obsowin <- as.owin(obspoly) #only step that requires maptools
unitname(obsowin) <- c("metre", "metres") #manually set units
plot(obsowin,
     main = "Observation Window as owin",axes=TRUE)


# # # # # # # # # # # # # # # # # # # # # # # # # # #
# 3. Extracting Raster Data for Observation Window  #
# # # # # # # # # # # # # # # # # # # # # # # # # # #
library("raster") # reads in remotely sensed raster data in a wide variety of formats.
#### First unzip the example raster data
#(raster data was compressed in a zip to save space)
rsdatafilepath <- system.file("extdata/demorsraster.zip",
                              package="lacunaritycovariance")
rsfiles <- unzip(rsdatafilepath,exdir=tempdir())

#### Open raster data file (this does not read the raster data)
xidataset<-raster(rsfiles[[2]])
# if the above doesn't load you may need to try opening rsfiles[[1]]
# (this may just be due to a filename extention convention)

#### Reads in the smallest rectangle possible around the observation window
xidataset <- crop(xidataset,extent(obspoly))
plot(xidataset, main="Raster Map around Observation Window")
plot(add=TRUE, obsbdry, lwd = 3)


# # # # # # # # # # # # # # # # # # # # # # # # # #
# 4. Convert Raster Data into spatstat im Object  #
# # # # # # # # # # # # # # # # # # # # # # # # # #
xiimage <- as.im.RasterLayer(xidataset) # uses package maptools
unitname(xiimage) <- c("metre","metres") # manually set units
# remove raster data outside observation window:
xiimage[setminus.owin(Frame(xiimage), obsowin)] <- NA
plot(xiimage, col = c("grey", "black"), main = "im Object Ready for lacunaritycovariance")

# # # #
# End #
# # # #
# The xiimage has values of 1, 0 and NA. It is ready be analysed
# using lacunaritycovariance.


# Thank you,
# Kassel

par(opa)
