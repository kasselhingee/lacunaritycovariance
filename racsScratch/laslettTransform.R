#calculating Laslett's transform

laslettTransform <- function(xi){
  xi <- as.im(xi)
  xi[is.na(as.matrix(xi))] <- 0 
  #gradient image for whole thing
  grad <- xi[,1:(xi$dim[2]-1)]-xi[,-1]
  
  #find lower tangent points. A grad value of 1 means exiting xi, and grad value of -1 means entering xi
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
        xi[rowIndex,en+1] <- 2 #the -1 is to take into
      } #label as tangent
    }
  }

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
  return(list(points=pointPat,image=xi))
}




