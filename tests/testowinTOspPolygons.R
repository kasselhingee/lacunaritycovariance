#testing spatial polygon conversion

library(spatstat, quietly = TRUE)
library(stationaryracsinference, quietly = TRUE)

#polygons that are all have simple boundary
data(polygontest, package="stationaryracsinference")
tess <- quadrats(maptools::as.owin.SpatialPolygons(polygontest))
tilelist <- tiles(tess)
id <- data.frame(id = 1:length(tilelist),stringsAsFactors = FALSE)
tileauxdata <- data.frame(t(simplify2array(strsplit(names(tilelist), "[ , ]"))),stringsAsFactors=FALSE)
tileauxdata <- tileauxdata[,c(3,6)]
names(tileauxdata) <- c("TileRow", "TileCol")
tileauxdata$TileRow <- factor(tileauxdata$TileRow,ordered=TRUE) 
tileauxdata$TileCol <- factor(tileauxdata$TileCol,ordered=TRUE)
auxdata <- cbind(id,tileauxdata)

spdf1 <- as.spdf.owin(tilelist, tileauxdata, idcolumn = NULL, 
                     proj4string = sp::CRS(sp::proj4string(polygontest)), match.ID=FALSE)
spdf1
                     
spdf2 <- sp::SpatialPolygonsDataFrame(as.SpatialPolygons.tess(tess), auxdata, match.ID = FALSE)
sp::proj4string(spdf2) <- sp::CRS(sp::proj4string(polygontest))
spdf2

#######################################################################
#introduce a non-simply boundary
polyowin <- maptools::as.owin.SpatialPolygons(polygontest)
polyowin <- setminus.owin(polyowin,owin(xrange=c(377850,377855),yrange=c(6527140,6527145)))
tess <- quadrats(polyowin)
tilelist <- tiles(tess)
id <- data.frame(id = 1:length(tilelist),stringsAsFactors = FALSE)
tileauxdata <- data.frame(t(simplify2array(strsplit(names(tilelist), "[ , ]"))),stringsAsFactors=FALSE)
tileauxdata <- tileauxdata[,c(3,6)]
names(tileauxdata) <- c("TileRow", "TileCol")
tileauxdata$TileRow <- factor(tileauxdata$TileRow,ordered=TRUE) 
tileauxdata$TileCol <- factor(tileauxdata$TileCol,ordered=TRUE)
auxdata <- cbind(id,tileauxdata)

spdf1 <- as.spdf.owin(tilelist, tileauxdata, 
                      proj4string = sp::CRS(sp::proj4string(polygontest)), match.ID=FALSE)
spdf1

spdf2 <- sp::SpatialPolygonsDataFrame(as.SpatialPolygons.tess(tess), auxdata, match.ID = FALSE)
sp::proj4string(spdf2) <- sp::CRS(sp::proj4string(polygontest))
spdf2


#my code misses the interior boundary, but Adrian's code has preserved it :)