#calculating Laslett's transform

library(spatstat)
data(heather)

xi <- heather$coarse
head(xi$m)
row <- xi$m[1,]
row <- as.numeric(row)
gradient <- row[-1]-row[1:(length(row)-1)] #difference between pixel at i and pixel at i-1

plot(gradient)
lines(row)

rowb <-as.numeric(xi$m[2,])
gradientb <- rowb[-1]-row[1:(xi$dim[2]-1)]

#trying EBImage
library(EBImage)
xi <- Image(as.numeric(as.matrix(heather$coarse)),heather$coarse$dim)
display(xi)
#but can't see a nice way to do this differencing thing


xi <- as.im(heather$coarse)

laslettTransform <- function(xi){
  xi[is.na(as.matrix(xi))] <- 0 
  #gradient image for whole thing
  Y <- xi[,-1]
  grad <- xi[,1:(xi$dim[2]-1)]-Y
  
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
        cat("tangent found at ",rowIndex,",",en,"\n")
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
  tangentPoints[,2] <- tangentPoints[,2]*xi$xstep+xi$xrange[1]
  tangentPoints[,1] <- tangentPoints[,1]*xi$ystep+xi$yrange[1]  
  
  #get a logical matrix - should be an easier way!
  xiLT <- matrix(data=NA,nrow=dim(xi)[1],ncol=dim(xi)[2])
  xiLT[!is.na(as.matrix(xi))] <- TRUE
  pointPat <- ppp(tangentPoints[,2],tangentPoints[,1],xrange=xi$xrange,yrange=xi$yrange,mask=xiLT)
  return(pointPat)
}

ltxi <- laslettTransform(xi)
layout(matrix(c(1,2),ncol=2,nrow=1))
plot(xi)
plot(ltxi)


