#given a polygon, function works out which UM tile is needed and returns raster data

require(raster)
#require(maptools)
require(rgdal)

#test rgdal's spTransform using two versions of poly03 in different CRS
test_reprojection <- function(){
  #polygon with same projection as UM data
  polyOGR_MGA50 <- readOGR("data","poly03") 
  coords_MGA50 <- polyOGR_MGA50@polygons[[1]]@Polygons[[1]]@coords  #an array of the coordinates of the polgon
  
  
  #same polygon with different projection 
  polyOGR_albers <- readOGR("data","poly03_albers") #in CRS GDA94 / australian albers
  #obsolete: readOGR loads polygon too :) polygon_albers <- readShapeSpatial("data/poly03_albers.shp",proj4string = crs(proj4string(polyOGR_albers)))
  coords_albers <- polyOGR_albers@polygons[[1]]@Polygons[[1]]@coords  #an array of the coordinates of the polgon
  polygon_albers_trans <- spTransform(polyOGR_albers,CRS("+proj=utm +zone=50 +south +ellps=GRS80 +units=m +no_defs")) #function from rgdal
  coords_albers_trans <- polygon_albers_trans@polygons[[1]]@Polygons[[1]]@coords #an array of the coordinates of the polgon
  
  #after transformation these should be the same
  return(abs(coords_albers_trans - coords_MGA50) < 1e-7)
}


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

