#' @title Functions for converting lists of im Ojects to on-disk rasterLayer objects
#' @rdname asondiskrasters
#' @export as.rasterlayer.im  offram imtoondiskras  imstoondiskras  rasterstoims
#' @description Functions for saving spatstat's \code{im} objects as \code{raster} objects out of RAM. Useful when computer has limited RAM and a fast solid state disk to store the images off RAM.
#' \code{as.rasterlayer.im} Converts a \code{spatstat} \code{im} into a \code{Raster} \code{rasterLayer}

#' @details These functions are useful when experiments generate a large collection of \code{im} objects.
#' The package \code{raster} has many useful tools for working with these images outside RAM.
#' I think it is wise to garbage collect \code{gc()} just after overwriting the im objects.

#' @param im A spatstat image object
#' @param layername Optional layer name to give to the created rasterLayer

#' @examples
#' xi <- heather$coarse
#' cvchats <- balancedracscovariances(xi, Frame(xi), modifications = "all")
#' cvchatsondisk <- imstoondiskras(cvchats)
#' cvchaatsbackasim <- rasterstoims(cvchatsondisk)

as.rasterlayer.im <- function(im, layername = NULL){
  stopifnot(class(im) == "im")
  imras <- raster::raster(im)
  if (!is.null(layername)){names(imras) <- layername}
  return(imras)
}

#' @param rasterobj A Raster* object
#' @param filenamesave The filename to save \code{rasterobj} to.
#' @describeIn asondiskrasters Create a rasterLayer (or other Raster* object) that is saved on disk out of another raster object.
#' If \code{filenamesave} is \code{NULL} then a temporary location is used.
offram <- function(rasterobj, filenamesave = NULL, overwrite = FALSE){
  if (is.null(filenamesave)){filenamesave <- tempfile()}
  offramrasterobj <- raster::writeRaster(rasterobj, filename = filenamesave, overwrite = overwrite)
  if (raster::fromDisk(offramrasterobj)){
    rm(rasterobj)
    gc()
    return(offramrasterobj)
  } else { stop("Raster object not saved off disk") }
}

#' @describeIn asondiskrasters Create a rasterLayer (or other Raster* object) that is saved on disk out of a \code{im} object.
imtoondiskras <- function(im, layername = NULL, filenamesave = NULL, overwrite = overwrite){
  imras <- as.rasterlayer.im(im, layername = layername)
  imras <- offram(imras, filenamesave = filenamesave, overwrite = overwrite)
  return(imras)
}

#' @param alist A list containing some \code{im} objects.
#' @param basefilename The filename that is used as root of all filenames. 
#' @param recursive If TRUE sublists of \code{alist} will also be modified.
#' @param overwrite If TRUE files, if they exist, will be overwritten.
#' @describeIn asondiskrasters Replaces all \code{im} elements in a list with on-disk \code{rasterLayer} objects.
#' Raster data is saved in \code{basefilename_[listentryname]}
imstoondiskras <- function(alist, basefilename = NULL, recursive = FALSE, overwrite = FALSE){
  class(alist) <- "list" #removes all other class attributes - for example stops alist checking if entries are ims
  if (is.null(basefilename)){basefilename <- tempfile()}
  noname <- vapply(names(alist), is.null, FALSE)
  if(sum(noname) > 0) {names(alist)[noname] <- 1:sum(noname)}
  if (recursive){
    arelist <- vapply(alist, function(x) "list" %in% class(x), FALSE) #ims should not have class list
    if (sum(arelist) > 0){
      alist[arelist] <- mapply(imstoondiskras, alist[arelist],
             basefilename = paste0(basefilename,"_",names(alist[arelist])),
             recursive = TRUE,
             overwrite = overwrite,
             SIMPLIFY = FALSE)
    }
  }
  areim <- vapply(alist, is.im, FALSE)
  if(sum(areim) == 0){return(alist)}
  alist[areim] <- 
    out <- mapply(imtoondiskras,
                         im = alist[areim],
                         layername = names(alist[areim]),
                         filenamesave = paste0(basefilename, "_", names(alist[areim])),
                         overwrite = overwrite,
                         SIMPLIFY = FALSE
                         )
  return(alist)
}

#' @describeIn asondiskrasters Replaces all code{rasterLayer} objects in \code{alist} with spatstat \code{im} objects. Uses MAPTOOLS
rasterstoims <- function(alist, recursive = FALSE){
  if (recursive){
    arelist <- vapply(alist, function(x) "list" %in% class(x), FALSE) #ims should not have class list
    if (sum(arelist) > 0){
      alist[arelist] <- mapply(rasterstoims, alist[arelist],
             recursive = TRUE,
             SIMPLIFY = FALSE)
    }
  }
  israsterlayer <- vapply(alist, function(x) "RasterLayer" %in% class(x), FALSE)
  if (sum(israsterlayer) == 0){return(alist)}
  alist[israsterlayer] <- lapply(alist[israsterlayer], maptools::as.im.RasterLayer)
  return(alist)
}
