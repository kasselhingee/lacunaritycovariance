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


xi <- as.im(heather$fine)

laslettTransform <- function(xi){
  xstep <- xi$xstep #saving this info for later
  ystep <- xi$ystep
  xrange <- xi$xrange
  yrange <- xi$yrange
  dimen <- xi$dim
  xiunits <- xi$units
  
  Y <- as.matrix(xi)
  Y[is.na(Y)] <- 0 
  #gradient image for whole thing
  grad <- Y[,1:(dim(xi)[2]-1)]-Y[,-1]

  #find lower tangent points. A grad value of 1 means exiting xi, and grad value of -1 means entering xi
  tangentPoints <- data.frame()
  for (rowIndex in 2:dim(xi)[1]){
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
        xi[rowIndex,en+1] <- 0 #because the proof in Cressie includes lower tangent points in the discretised background
        tangentPoints <- rbind(tangentPoints,c(rowIndex,en+1))
        } #label as tangent
    }
  }
  colnames(tangentPoints) <- c("x","y")
  
  #remove the stuff in xi, keeping the tangent points - because thats the discretisation used in the proof in Cressie's book
  for (rowIndex in 1:dim(xi)[1]){
    rowVal <- xi[rowIndex,]
    rowValinxi <- (rowVal==1)
    rowValNew <- c(rowVal[!rowValinxi],rep(NA,sum(rowValinxi)))
    xi[rowIndex,] <- rowValNew
  }
  xiLT <- owin(mask = as.logical(xi+1), unitname=xiunits)
  
  return(c(xi)
}

ltxi <- laslettTransform(xi)
layout(matrix(c(1,2),ncol=2,nrow=1))
plot(xi)
plot(ltxi)


