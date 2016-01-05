library(maptools) #with rgeos installed (for ubuntu need to install GEOS separately - waiting to see if needed!)
library(raster)
library(spatstat)

#read boundary in and convert to spatstat window
boundaryMAPTOOLS <- readShapeSpatial("data/poly01_polygons.shp")
boundaryMAPTOOLS <- readShapeSpatial("../../Data/CBC_Feed_areas_requiring_investigation_split.shp")
boundaryOWIN <- as.owin(boundaryMAPTOOLS[16408,1]) #this is actually row 16407 in the QGIS viewer
boundaryOWIN <- as.owin(boundaryMAPTOOLS)  #

##read in class map convert to spatstat window
XiRASTER <- raster("data/inputBinaryMap.ers")
XiRASTER <- raster("C:/CCI-02_Work/processing_102/UM2009/2035SW_moorSE/2009_mar_UM_ucd03_2035SW_moorSE_gda94_mga50_tre_raw_CSIRO_2012092313_ver01.ers")
XiRASTER <- crop(XiRASTER,boundaryMAPTOOLS[16408,1]) #this process takes some time, but doesn't actually create anything in memory
plot(XiRASTER)
plot(boundaryOWIN,add=TRUE)
XiIMAGE <- as.im(XiRASTER) #the data is now saved into memory. It is pretty inefficient - 149MB for a 39MB ERMapper file
XiOWIN <- as.owin(XiIMAGE)

#check that area of a mask owin is number of pixels in class
area.owin(XiOWIN)
sum(XiOWIN$m)*XiOWIN$xstep*XiOWIN$ystep  #also works to use as.matrix(XiOWIN)


#estimate exact variance
varPexactEst(covprobEstimate,covarianceMap,boundaryOWIN)


maxXshiftdistance <- 3.1
c(-rev(seq(0,maxXshiftdistance,by=0.2)),seq(0.2,maxXshiftdistance,by=0.2))



plot(shift.owin(Xiinside,vec=shiftVector))
plot(boundaryOWIN,add=TRUE)
plot(rData)
plot(polyWindowowin,add=TRUE)


###############heather data
data(heather)
XiOWIN <- heather$coarse
windowOWIN <- owin(xrange=heather$coarse$xrange,yrange=heather$coarse$yrange)

coverageProb <- covpest(XiOWIN,windowOWIN)

covariancePt <- covarianceEstAtPoint(XiOWIN,windowOWIN,c(3,5))

covarianceDirectEst <- covarianceMapEst_direct(XiOWIN,windowOWIN,1,1)
filled.contour(covarianceDirectEst$covariance)


covarianceFcn <- covarianceRACS(XiOWIN,windowOWIN)
plot(covarianceFcn$covariance)
plot(covarianceFcn$covariance - coverageProb*coverageProb)

Z <- HestInPolygonalWindow(XiOWIN,windowOWIN)
plot(Z)
plot(add=TRUE,Hest(XiOWIN))#since window is rectangular should be the same


###########################
#simulating a boolean model of discs of random radius with lognormal distribution
w <- owin(xrange=c(0,10),yrange=c(0,10))
lambda <- 1
meanlog <- -1
sdlog <- 0.5
plot(0.01*(1:200),plnorm(0.01*(1:200),meanlog=meanlog,sdlog=sdlog),type="l")
#probability of getting a radius above 2 is less than 0.0004
r <- 2 #a dilation distance chosen such that the probablity of a germ with a grain that overlaps w is very small


wsim <- dilation(w,r)
#have to simulate in a much larger area than the observation window (because grains with centres outside the window should still be observed)
pp <- rpoispp(lambda,win=wsim,nsim=1,drop=TRUE)
plot(pp)
#need to make this work on multiple simulations? so I can make many all at once?

radius <- rlnorm(pp$n,meanlog=meanlog,sdlog=sdlog) #prepare a random radius for each point

pointlocations <- cbind(X=pp$x,Y=pp$y)
pointlocations <- split(cbind(pointlocations),row(pointlocations)) #split matrix into a list of the rows
grains <- mapply(disc,radius = radius,centre=pointlocations,SIMPLIFY=FALSE) #calculate grains with their locations

#take union of all grains
xisim <- union.owin(as.solist(grains))
plot(wsim)
plot(add=TRUE,xisim)
plot(add=TRUE,w)

#intersect back to get the observation window
xi <- intersect.owin(xisim,w)
plot(w)
plot(add=TRUE,xi)

xilist <- solist(xi,
                 rboollognormdiscs(w,2,1,-1,0.5),
                 rboollognormdiscs(w,2,1,-1,0.5),
                 rboollognormdiscs(w,2,1,-1,0.5),
                 rboollognormdiscs(w,2,1,-1,0.5),
                 rboollognormdiscs(w,2,1,-1,0.5)
                 )
plot(xilist,nrow=2)


#function for attaching grains from a library of grains!
grainlib <- as.solist(mapply(disc,radius = rlnorm(5,meanlog=-1,sdlog=0.5),SIMPLIFY=FALSE) )

#dummy pp
pp <- rpoispp(0.05,win=wsim,nsim=1,drop=TRUE)

#shift and attach grains randomly
grains <- sample(grainlib,size=pp$n,replace=TRUE) 
pointlocations <- cbind(X=pp$x,Y=pp$y)
pointlocations <- split(cbind(pointlocations),row(pointlocations)) #split matrix into a list of the rows
shiftedgrains <- mapply(shift.owin,grains,vec=pointlocations,SIMPLIFY=FALSE)
uniongrains <- union.owin(as.solist(shiftedgrains))

plot(wsim)
plot(uniongrains,add=TRUE)
plot(pp,add=TRUE,pch="+")

xi <- placegrainsfromlib(pp,grainlib)
plot(wsim)
plot(xi,add=TRUE)

plot(grainlib,nrow=2,equalscales=TRUE)



#spectral density example, use boolean model of discs
#Generating a germ-grain models where germs are a Poisson Point process, and grains are 2 or 3 different disc sizes.
grainlib <- solist(disc(radius=10))
bufferdist <- 12 #chosen to be larger than the largest radius in library
w <- owin(xrange=c(0,500),yrange=c(0,500)) #large numbers of points makes attaching grains take forever!?
pp <- rpoispp(lambda=2.2064E-3,win=dilation(w,bufferdist),nsim=1,drop=TRUE)
xibuffer <- placegrainsfromlib(pp,grainlib)
xi <- intersect.owin(xibuffer,w)
plot(xi)
#40seconds with w <- owin(xrange=c(0,500),yrange=c(0,500)) and lambda=2.2064E-3, 
#using mask values seems 2+ times slower - would need to write a whole lot of code to discretise the whole process to make it faster"

p <- covpest(xi,Frame(xi))
X <- as.im(xi,dimyx=c(511,511))


specdens <- spectraldensity(xi,w,dimyx=c(512,512))
plot(specdens,clipwin=owin(xrange=c(-50,50),yrange=c(-50,50)))
plot(specdensim,clipwin=owin(xrange=c(-50,50),yrange=c(-50,50)))

#compare to covariance function (keep this as a test later)
covar <- covarianceRACS(xi,Frame(xi))
plot(covar$covariance-p*p)#this is kinda what I'd expect fairly uniform. Although its got some rectangular effect happening
plot(covar$covariance-p*p,clipwin=owin(xrange=c(-20,20),yrange=c(-20,20))) #a spike kinda like a radius of 10 as discs are all size 10
M <- covar$covariance$v
M[is.na(M)] <- 0 
specdensB <- fft(M)/length(M)

image(Re(specdensB[110:140,110:140]))
#doesn't fit!!

#test on a boolean model of rectangles
areaGrainRect = 314.16
sideRatio = 2   #area is sideRatio*sideAlength*sideAlength
sideAlength = sqrt(areaGrainRect/sideRatio)
grainlib <- solist(owin(xrange = c(0,sideAlength),yrange=c(0,sideAlength*sideRatio)))
bufferdist <- 2*sideAlength*sideRatio #chosen to be larger than the largest radius in library
w <- owin(xrange=c(0,500),yrange=c(0,500)) #large numbers of points makes attaching grains take forever!?
pp <- rpoispp(lambda=2.2064E-3,win=dilation(w,bufferdist),nsim=1,drop=TRUE)
xibuffer <- placegrainsfromlib(pp,grainlib)
xi <- intersect.owin(xibuffer,w)
plot(w)
plot(xi,add=TRUE)
X <- as.im(xi,dimyx=c(512,512))
covpest(xi,w)
specden <- spectraldensity(xi,w,dimyx=c(512,512))
plot(specden,clipwin=owin(xrange=c(-50,50),yrange=c(-50,50)))

#continuing onto kernel smoothing anyway. package fields does exactly what I want using the image() function!
#it probably isn't too hard to write my own using FFT
library(fields)
EpanechnikovFcnFields <- function(sz){#WARNING: operates on a vector
  if (sz>1) {return(0)}
  if (sz <= 1) {return(2/pi*(1-sz))}
  return(NA)
}
EpanechnikovFcn <- function(X,Y){#WARNING: operates in 2D only on a vector of things 
  stopifnot(length(X)==length(Y))
  result <- vector(length=length(X),mode="numeric")
  sz <- sqrt((X*X)+(Y*Y))
  result[sz>1] <- 0
  result[sz<=1] <- 2/pi*(1-sz[sz<=1]^2)
  return(result)
}
bandwidth = 6
xstep = specdens$xstep
ystep = specdens$ystep
supportwidth = bandwidth
X <- seq(0,supportwidth*1.5+xstep,by=xstep) #much larger than support width to avoid boundary issues?
Y <- seq(0,supportwidth*1.5+ystep,by=ystep)
mat <- outer(X/bandwidth,Y/bandwidth,FUN="EpanechnikovFcn") #rows correspond to xstep - just a quirk of outer!
mat <- t(mat) #columns correspond to changes in X, rows correspond to changes in Y!
#reflect out to all corners
mat <- mat[,c((ncol(mat)):2,1:ncol(mat))]
mat <- mat[c((nrow(mat)):2,1:nrow(mat)),]
kernelfcn <- im(mat,xcol=c(-X[length(X):2],X),yrow=c(-Y[length(Y):2],Y))
#apply convolve.im
smspecdens <- convolve.im(specdens,kernelfcn)
plot(smspecdens,clipwin=owin(xrange=c(-50,50),yrange=c(-50,50)))
plot(smspecdens,clipwin=owin(xrange=c(-100,100),yrange=c(-100,100)))

smspecdens <- spectraldensity(xi,Frame(xi),6,dimyx=c(512,512))
plot(smspecdens)
dim(smspecdens)
plot(smspecdens$xcol, smspecdens[257,])
#covariance estimate
ifspecdens <- fft(as.matrix(smspecdens),inverse = TRUE)/((2*pi)^2)
image(Im(ifspecdens))
image(Re(ifspecdens)[0:50,0:50])
image(Re(ifspecdens)[0:50,0:50])
plot(Re(ifspecdens)[1,1:20])

#theoretical covariance of a Boolean model is given (with proof) in Eqn 3.18 of Chiu 2014. Examples for particular grains can be found in B"ohm 2002
thcovDeterministicDiscs <- function(r,lambda,discr){
  if (r>=2*discr){expectedsetcovariance <- 0}
  else {
    expectedsetcovariance <- 2*discr^2*acos(r/(2*discr)) - (r/2)*sqrt(4*discr^2-r^2)
  }
  p <- 1-exp(-pi*discr^2*lambda)
  covariance <- 2*p-1+(1-p)^2*exp(lambda*expectedsetcovariance)
}
plot(1:100/5,lapply(1:100/5,thcovDeterministicDiscs,lambda = lambda, discr=discr))

thcovDeterministicDiscs_vec <- function(X,Y,lambda,discr){
  rlist <- sqrt(X^2+Y^2)
  covar <- vector(length(rlist),mode="numeric")
  for (i in 1:length(rlist)){
    covar[i] <- thcovDeterministicDiscs(rlist[i],lambda=lambda,discr=discr)
  }
  return(covar)
}
plot(1:100,thcovDeterministicDiscs_vec(1:100,1:100,lambda=lambda,discr=discr))
mat <- outer(0:250/10,0:250/10,FUN="thcovDeterministicDiscs_vec",lambda=lambda,discr=discr) #rows correspond to xstep - just a quirk of outer!
image(mat)
mat <- t(mat) #columns correspond to changes in X, rows correspond to changes in Y!
#reflect out to all corners
mat <- mat[,c((ncol(mat)):2,1:ncol(mat))]
mat <- mat[c((nrow(mat)):2,1:nrow(mat)),]
image(mat)

thspecdens <- fft(mat)/length(mat)
dim(thspecdens)
plot(1:100,Re(thspecdens[1,1:100]))
plot(1:100,Im(thspecdens[1,1:100]))
nr = nrow(thspecdens)
nc = ncol(thspecdens)
thspecdens <- thspecdens[ ((-(nr-1)/2):((nr-1)/2)) %% (nr) + 1,((-nc/2):(nc/2)) %% (nc) + 1]
filled.contour(Re(thspecdens))

plot(specdens,clipwin=owin(xrange=c(-50,50),yrange=c(-50,50)))

smspecdensFields <- image.smooth(as.matrix(specdens),
                           kernel.function=EpanechnikovFcnFields,
                           dx=specdens$xstep, #passing dimensions
                           dy=specdens$ystep, #passing dimensions
                           xwidth = 6, #zero padding around outside
                           theta=6)  #this is the bandwidth
image.plot(smspecdens)
plot(as.im(smspecdens,W=specdens))


#test on a known function. (A cosine wave)
xcol <- 0:15/16
arow <- cos(3*xcol*2*pi)
plot(arow)
Mim <- as.im(matrix(arow,nrow=16,ncol=16,byrow=TRUE),xcol=xcol,yrow=xcol)
plot(Re(fft(arow)/16))
#has value 0.5 at 4 and 14. Is this correct?
#fourier transform evaluated at 3 Hz should be the only spike. 
#the value at location 4 corresponds to a frequency of 3/16 (periodic boundary effect?)
#the value at location 14 corresponds to a frequency of 13/16
specdenT <- (Re(fft(arow))^2+Im(fft(arow))^2)/((16*16)^3)
specdenT <- specdenT[ ((-16/2):(16/2)) %% (16) + 1]
specdens <- spectraldensity(Mim)
plot(specdens)
specdens
##I still not sure if is correct. What are the two spikes for!??

#apply to heather
xi <- heather$coarse
unsmspecdens <- unsmoothedspectraldensity(xi,Frame(xi))
plot(unsmspecdens)
specdens <- spectraldensity(xi,Frame(xi),20)
plot(specdens)

plot(solist(smoothed = specdens,unsmoothed = unsmspecdens),
     main = "Spectral Density Estimates",
     equal.scales=TRUE)

#for isotropic RACS can take average over angles



#testing out varblock on Hest
varblock(heather$coarse,fun=Hest)
#only works on ppp data currently! but there are lots of tools that (eg quadrats() ) that varblock uses and will be useful to me too
