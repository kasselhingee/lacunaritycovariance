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


#try out spatstats convolve (FFT)
xi <- as.im(heather$coarse)
xi[is.na(as.matrix(xi))] <- 0 #possibly unsupported use of spatstat
#gradient image for whole thing
Y <- xi[,-1]
grad <- xi[,1:(xi$dim[2]-1)]-Y
plot(grad)
plot(xi)

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
      xi[rowIndex,en] <- 2 #the -1 is to take into
      } #label as tangent
  }
}

plot(xi)
