#calculating Laslett's transform

#' @title Laslett's Transform
#' @export lasletttransform findlowlefttangentpts
#' 
#' @description Performs Laslett's tranform on a raster map. \code{findlowlefttangentpts} is function used to find the tangent points.
#' See ** reference.
#' 
#' @param xi An owin mask object.
#' @return \code{lasletttranform} Returns a ppp of the tranformed location of tangent points with transformed window.
#'        \code{findlowlefttangentpts} returns two lists, labelled X and Y respectively of the X and Y coordinates of tangent point.

#' @examples 
#' 
#' xi <- rotate.owin(heather$coarse)
#' ltxi <- lasletttransform(xi)
#' plot(ltxi)
#' 
#' 
#' #test on a true boolean model
#' xi <- rBooleanDetermDiscs(2.2064E-3,10,owin(xrange=c(0,500),yrange=c(0,500)))
#' xim <- as.mask(xi,eps=c(1,1))
#' ltxi <- lasletttransform(xi)
#' Ks <- Kest(ltxi)
#' plot(Ks)
#' 
lasletttransform <- function(xi){
  tangentpointslist <- findlowlefttangentpts(xi)
  
  xi <- as.im(xi)
  xi[is.na(as.matrix(xi))] <- 0 
  #use the pixels containing these points as the tangent points (will it actually matter if this is slightly off - it all approximates the same thing as resolution gets infintely fine)
  xi[ppp(tangentpointslist$X,tangentpointslist$Y,window=Frame(xi))] <- 2
  
  #remove the stuff in xi, keeping the tangent points - because thats the discretisation used in the proof in Cressie's book
  for (rowIndex in 1:xi$dim[1]){
    rowVal <- xi[rowIndex,]
    rowValinxi <- (rowVal==1)
    rowValNew <- c(rowVal[!rowValinxi],rep(NA,sum(rowValinxi)))
    xi[rowIndex,] <- rowValNew
  }
  #find tangent points in converted image
  tangentPoints <- arrayInd(which(as.matrix(xi)==2),dim(xi))
  #convert tangent points to cartesian coordinates
  tangentPoints[,2] <- (tangentPoints[,2]-1)*xi$xstep+xi$xcol[1] #to get the correct coordinates!
  tangentPoints[,1] <- (tangentPoints[,1]-1)*xi$ystep+xi$yrow[1]  
  
  #get a logical matrix - should be an easier way!
  xiLT <- matrix(data=NA,nrow=dim(xi)[1],ncol=dim(xi)[2])
  xiLT[!is.na(as.matrix(xi))] <- TRUE
  xiLT[is.na(as.matrix(xi))] <- FALSE
  pointPat <- ppp(tangentPoints[,2],tangentPoints[,1],xrange=xi$xrange,yrange=xi$yrange,mask=xiLT)
  return(points=pointPat)
}
#' @describeIn lasletttransform Returns a list of X and Y coordinates of tangent points
#' 
findlowlefttangentpts <- function(xi){
    xi <- as.im(xi)
    xi[is.na(as.matrix(xi))] <- 0 
    #gradient image for whole thing
    grad <- xi[,1:(xi$dim[2]-1)]-xi[,-1]
    
    #find lower tangent points. A grad value of 1 means exiting xi, and grad value of -1 means entering xi
    tangentpoints=list(x = vector(),y = vector())
    for (rowIndex in 2:xi$dim[1]){
      rowa <- grad[rowIndex,]#rows start at the bottom!!
      rowb <- grad[rowIndex-1,]
      rowaexits <- which(rowa==1)
      rowaenter <- which(rowa==-1)
      
      rowbexits <- which(rowb==1)
      rowbenter <- which(rowb==-1)
      rowbchanges <- which((rowb==1) | (rowb==-1))
      #search for an entrance,exit in rowa that has no entrance or exit in rowb over the same length
      for (en in rowaenter){
        ex <- rowaexits[rowaexits>en][1]
        if (length(rowaexits[rowaexits>en])==0){next}
        tangent <- ((min((rowbenter< en) | (rowbenter >= ex))) && #if no enter point along interval (this is 4 nbhd because of strictly less/greater than) (#something better than max??)
                      (xi[rowIndex-1,en] == 0))    #and first pixel isn't in set 
        if(tangent){
          tangentpoints$x = append(tangentpoints$x,en+1)
          tangentpoints$y = append(tangentpoints$y,rowIndex)
        } 
      }
    }
    #convert list of tangent points to spatial locations
    return(list(X = ((tangentpoints$x-1) * xi$xstep)+xi$xcol[1],
                Y = (xi$ystep * (tangentpoints$y-1))+xi$yrow[1]))
}
