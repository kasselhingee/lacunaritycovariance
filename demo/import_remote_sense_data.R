
opa <- par(mfrow=c(1,1))

##################
# IMPORTING REMOTE SENSING DATA DEMO
# By Kassel Liam Hingee
##################

# This document demonstrates how to convert raster data stored in a
# remote sensing format and an observation window in an ESRI shapefile
# format into a format suitable for the package lacunaritycovariance.

# For this demo you will need the following R packages installed: 
#  raster (and terra)
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

# # # # # # # # # # # # # # # # # # # # # #
# 1. Reading ESRI Shapefile in SF Object  #
# # # # # # # # # # # # # # # # # # # # # #
regionfilepath <- system.file("extdata/aregionofinterest.shp", package="lacunaritycovariance")
obspoly <- sf::st_read(regionfilepath)
#print the coordinate projection of the polygon data for sanity
plot(obspoly, main = "Observation Window as a SF Object")


# # # # # # # # # # # # # # # # # # # # # # # # # # #
# 2.Converting Observation Window into owin Format  #
# # # # # # # # # # # # # # # # # # # # # # # # # # #
obsowin <- as.owin(obspoly) 
unitname(obsowin) <- c("metre", "metres") #manually set units
plot(obsowin,
     main = "Observation Window as owin",axes=TRUE)


# # # # # # # # # # # # # # # # # # # # # # # # # # #
# 3. Extracting Raster Data for Observation Window  #
# # # # # # # # # # # # # # # # # # # # # # # # # # #
library("terra") # reads in remotely sensed raster data in a wide variety of formats.
#### First unzip the example raster data
#(raster data was compressed in a zip to save space)
rsdatafilepath <- system.file("extdata/demorsraster.zip",
                              package="lacunaritycovariance")
rsfiles <- unzip(rsdatafilepath,exdir=tempdir())

#### Open raster data file (this does not read the raster data)
xidataset<-terra::rast(rsfiles[[2]])
# if the above doesn't load you may need to try opening rsfiles[[1]]
# (this may just be due to a filename extention convention)

#### Reads in the smallest rectangle possible around the observation window
xidataset <- crop(xidataset,obspoly)
plot(xidataset, main="Raster Map around Observation Window")
plot(add=TRUE, obspoly, lwd = 3, col = NA)


# # # # # # # # # # # # # # # # # # # # # # # # # #
# 4. Convert Raster Data into spatstat im Object  #
# # # # # # # # # # # # # # # # # # # # # # # # # #
# As of Oct 2023, the below is the best method I'm aware of for converting to a spatstat image object (the previous method was via the now obsolete maptools package).
xiimage <- spatstat.geom::im(
  as.matrix(xidataset, wide = TRUE)[nrow(xidataset):1, ],
   xrange = sf::st_bbox(xidataset)[c("xmin", "xmax")],
   yrange = sf::st_bbox(xidataset)[c("ymin", "ymax")])
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
