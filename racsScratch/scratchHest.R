##loading a shape file into sp format
#requires rgeos?, maptools, spatstat, raster, rgdal,

library(maptools) #with rgeos installed (for ubuntu need to install GEOS separately - waiting to see if needed!)

polyWindow <- readShapeSpatial("data/poly01_polygons.shp")
y <- readShapeSpatial("data/rect01_polygons.shp")
class(polyWindow)
plot(polyWindow)
class(y)
plot(y)

##convert to spatstat polygonal region
library(spatstat)
polyWindowowin <- as.owin(polyWindow)
plot(polyWindowowin)


##read in raster data within polygon
library(raster) #with rgdal installed
rData <- raster("data/inputBinaryMap.ers")
rData <- as.im(rData)
rasterD <- as.owin(rData)
plot(rData)
plot(polyWindowowin,add=TRUE)



RasterInsidePolyWindow <- intersect.owin(rasterD,polyWindowowin)
plot(RasterInsidePolyWindow)



##Alter Hest!? Copied from Hest code in spatstat, trying to replicate. Using code on github!
#The below works perfectly when using rasterD, and the range of rasterD as the boundary: the stuff below gives identical results to Hest of rasterD. 
condition <- function(x) { (x - x[1])/(1-x[1]) }
reconstitute <- function(x, p) { p + (1-p) * x }


X <- rasterD
#boundaryIn <- owin(rasterD$xrange,rasterD$yrange)
boundaryIn <- polyWindowowin
#X <- rasterD



r = NULL
breaks = NULL
correction = c("km", "rs", "han")
#correction <- c("rs", "km", "cs")  #BUG ORIGINAL - should be km, rs and han
conditional = TRUE

if(!(is.ppp(X) || is.psp(X) || is.owin(X)))
  stop("X should be an object of class ppp, psp or owin")
## handle corrections
if(is.null(correction))
  correction <- c("rs", "km", "cs")
correction <- pickoption("correction", correction,
                         c(none="none",
                           raw="none",
                           border="rs",
                           rs="rs",
                           KM="km",
                           km="km",
                           Kaplan="km",
                           han="han",
                           Hanisch="han",
                           best="km"),
                         multi=TRUE)
corxtable <- c("km", "rs", "han", "none") 
corx <- as.list(corxtable %in% correction)
names(corx) <- corxtable
## compute distance map
D <- distmap(X)
#B <- attr(D, "bdry") #for X this is a distance map from the rectangular boundary
B <- distmap.owin(boundaryIn,dimyx=X$dim,invert=TRUE) #warning! no gaurantee the pixel locations are the same as in the raster object?
#B <- distmap.owin(polyWindowowin,dimyx=X$dim,invert=TRUE) #warning! no gaurantee the pixel locations are the same as in the raster object?
# X$dim is in ny, nx order
W <- as.owin(D)
## histogram breakpoints 
dmax <- summary(D)$max
breaks <- handle.r.b.args(r, breaks, W, NULL, rmaxdefault=dmax)
rval <- breaks$r
##  extract distances and censoring distances
dist <- as.vector(as.matrix(D))
bdry <- as.vector(as.matrix(B))
ok <- !is.na(dist) && !is.na(bdry)
dist <- dist[ok]
bdry <- bdry[ok]
## delete zero distances #KH and distances for pixels outside the window
if(is.owin(X)) {
  pos <- (dist > 0) #pixels outside raster Xi
  inWindow <- ( bdry>=0 ) #KH add, includes pixels on edge of window, otherwise mean(pos) is different for rectangular regions
  #areafraction <- 1 - mean(pos)
  areafraction <- 1 - (sum(pos & inWindow)/sum(inWindow))
  #dist <- dist[pos]
  dist <- dist[pos & inWindow] #KH
  #bdry <- bdry[pos]
  bdry <- bdry[pos & inWindow] #KH
}
## censoring indicators
d <- (dist <= bdry)
##  observed distances
o <- pmin.int(dist, bdry)
## calculate estimates
Z <- censtimeCDFest(o, bdry, d, breaks,
                    KM=corx$km,
                    RS=corx$rs,
                    HAN=corx$han,
                    RAW=corx$none,
                    han.denom=if(corx$han) eroded.areas(W, rval) else NULL,
                    tt=dist)
## conditional on d > 0 ?
if(is.owin(X)) {
  if(conditional) {
    if(corx$km)   Z$km  <- condition(Z$km)
    if(corx$rs)   Z$rs  <- condition(Z$rs)
    if(corx$han)  Z$han <- condition(Z$han)
    if(corx$none) Z$raw <- condition(Z$raw)
  } else {
    if(corx$km)   Z$km  <- reconstitute(Z$km, areafraction) 
    if(corx$rs)   Z$rs  <- reconstitute(Z$rs, areafraction) 
    if(corx$han)  Z$han <- reconstitute(Z$han, areafraction) 
    if(corx$none) Z$raw <- reconstitute(Z$raw, areafraction) 
  }
}
## relabel
Z <- rebadge.fv(Z, substitute(H(r), NULL), "H")
unitname(Z) <- unitname(X)

plot(Z)
Ht <- Hest(rasterD)
plot(Ht,add=TRUE)

source("HestInPolygonalWindow.R")
Z2 <- HestInPolygonalWindow(RasterInsidePolyWindow,polyWindowowin)
Z <- HestInPolygonalWindow(rasterD,polyWindowowin)
plot(Z2)
plot(Z2-Z)
