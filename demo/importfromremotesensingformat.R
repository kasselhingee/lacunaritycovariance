


library("stationaryracsinference")
library("rgdal")
library("maptools")
library("raster")




regionfilepath <- system.file("extdata", package="stationaryracsinference")
obspoly <- readOGR(regionfilepath,"aregionofinterest")
#print the coordinate projection of the polygon data for sanity
crs(proj4string(obspoly)) 




obsbdry <- as.owin(obspoly) #only step that requires maptools
unitname(obsbdry) <- c("metre", "metres")
plot(obsbdry,
     main = "A Region of Interest / Observation Window",axes=TRUE)




#First unzip raster data
#(raster data was compressed in a zip to save space)
rsdatafilepath <- system.file("extdata/demorsraster.zip",
                              package="stationaryracsinference")
rsfiles <- unzip(rsdatafilepath,exdir=tempdir())
xidataset<-raster(rsfiles[[2]])
# if this doesn't load you may need to try opening rsfiles[[1]]
# (it may just be due to a filename extention convention)

xidataset <- crop(xidataset,extent(obspoly))
plot(xidataset,main="Tree Canopy Map and Region of Interest")
plot(add=TRUE,obsbdry)




xiimage <- as.im(xidataset)
unitname(xiimage) <- c("metre","metres")
xiimage[setminus.owin(Frame(xiimage), obsbdry)] <- NA
par(bg = "grey")
plot(xiimage, col = c("white", "black"))
