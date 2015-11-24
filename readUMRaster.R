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
#polygon with same projection as UM data
polyOGR_MGA50 <- readOGR("data","poly03") 
polygon_MGA50 <- readShapeSpatial("data/poly03.shp",proj4string = crs(proj4string(polyOGR_MGA50)))
coords_MGA50 <- polygon_MGA50@polygons[[1]]@Polygons[[1]]@coords  #an array of the coordinates of the polgon


#same polygon with different projection 
polyOGR_albers <- readOGR("data","poly03_albers") #in CRS GDA94 / australian albers
polygon_albers <- readShapeSpatial("data/poly03_albers.shp",proj4string = crs(proj4string(polyOGR_albers)))
coords_albers <- polygon_albers@polygons[[1]]@Polygons[[1]]@coords  #an array of the coordinates of the polgon
polygon_albers_trans <- spTransform(polygon_albers,CRS("+proj=utm +zone=50 +south +ellps=GRS80 +units=m +no_defs"))
coords_albers_trans <- polygon_albers_trans@polygons[[1]]@Polygons[[1]]@coords #an array of the coordinates of the polgon

#after transformation these should be the same
abs(coords_albers_trans - coords_MGA50) < 1e-7
  
  
  
polygon_albers == polygon_MGA50
plot(polygon_MGA50,col ="red")

proj4string(polygonMAPTOOLS) == "+proj=utm +zone=50 +south +ellps=GRS80 +units=m +no_defs"



#actual function
readUMraster <- function(polygon,prodCode,repository){
  if (is.na(proj4string(polygon))){ warning("Can not read pojection information of polygon - differences between projections of UM data and polygon will lead to importing the incorrect region")}
  else {
    if (proj4string(polygon) != "+proj=utm +zone=50 +south +ellps=GRS80 +units=m +no_defs") {stop("coordinate reference system of vector does not match UM reference system and no reprojection attempted")}
  }
  tileDims <- read.table("tileDimensions.txt",header=TRUE,stringsAsFactors=FALSE)
  extents <- extent(polygon)
  #find any tiles that completely contain extents rectangle
  tile <- tileDims[((tileDims$TLEASTING < extents[1]) 
                    & (tileDims$TLNORTHING > extents[4]) 
                    & (tileDims$BREASTING > extents[2]) 
                    & (tileDims$BRNORTHING < extents[3])),]
  if (dim(tile)[1] == 0 ){stop("no UM tile covered the polygon completely")}
  fileList <- list.files(paste(repository,tile$mapID[1],"/",sep=""))
  
  filename <- findFile(fileList,prodCode)[1]
  XiRASTER <- raster(paste(repository,"/",tile$mapID[1],"/",filename,".ers",sep=""))
  
  XiRASTER <- crop(XiRASTER,polygon)
}


#given a file list finds the ones corresponding to a particular product
#returns
findFile <- function(fileList,prodCode,procCode="any"){#"any" will return the most recent product irrespective of its processing
  fileName="ERROR"
  allDataFiles =fileList[!((grepl(".",fixed=TRUE,fileList))|grepl("patches",fileList)|grepl("old_versions",fileList))]
  if (length(allDataFiles)>1){
    splitDataNames=do.call("rbind",strsplit(allDataFiles,"_"))
    splitDataNames=data.frame(splitDataNames,stringsAsFactors=FALSE)
    splitDataNames[,10]=factor(splitDataNames[,10],c("9-7","non","nod","raw","cal","edt"),ordered=TRUE)
    splitDataNames[,12]=as.numeric(splitDataNames[,12])
    if (procCode!="any"){filesBool=((splitDataNames[,9]==prodCode)&(splitDataNames[,10]==procCode))
    } else {filesBool=(splitDataNames[,9]==prodCode)}
    if (max(filesBool)==1){#aka at least one appropriate file exists
      ordering=order(splitDataNames[filesBool,12],splitDataNames[filesBool,10],decreasing=TRUE)
      #files_index = which(splitDataNames[filesBool,12]==max(splitDataNames[filesBool,12])) # get most recent
      fileName = allDataFiles[filesBool][ordering[1]]
      return(fileName)
    }
  }
  
  if (length(allDataFiles)==1){
    splitDataNames=strsplit(allDataFiles[1],"_")
    if ((splitDataNames[[1]][9]==prodCode)&&((procCode=="any") || (splitDataNames[[1]][10]==procCode))){
      fileName = allDataFiles[1]
      return(fileName)
    }
  }
  
  cat(paste("Failed to find", prodCode,procCode,"in:",sep=" "),fileList,sep="\n")
  return("ERROR")	
}

