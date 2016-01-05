#testing spectral density calculations
library(stationaryracsinference)

#use boolean model of discs
#Generating a germ-grain models where germs are a Poisson Point process, and grains are 2 or 3 different disc sizes.
grainlib <- solist(disc(radius=10))
bufferdist <- 12 #chosen to be larger than the largest radius in library
w <- owin(xrange=c(0,500),yrange=c(0,500)) #large numbers of points makes attaching grains take forever!?
pp <- rpoispp(lambda=2.2064E-3,win=dilation(w,bufferdist),nsim=1,drop=TRUE)#lambda from B\"{o}m (2002) - chosen to make coverage probability very close to 0.5
xibuffer <- placegrainsfromlib(pp,grainlib)
xi <- intersect.owin(xibuffer,w)
plot(xi)
#40seconds with w <- owin(xrange=c(0,500),yrange=c(0,500)) and lambda=2.2064E-3, 
#using mask values seems 2+ times slower - would need to write a whole lot of code to discretise the whole process to make it faster


smoothedspecdens <- spectraldensity(xi,Frame(xi),6,dimyx=c(512,512))
plot(smoothedspecdens,clipwin=owin(xrange=c(-50,50),yrange=c(-50,50)))
plot(smoothedspecdens,clipwin=owin(xrange=c(-100,100),yrange=c(-100,100)))
plot(smspecdens$xcol, smspecdens[257,])

#theoretical covariance of a Boolean model is given (with proof) in Eqn 3.18 of Chiu 2014. Examples for particular grains can be found in B"ohm 2002
#set covariance of the grains:
setcovdisc <- function(r,discr){
  if (r>=2*discr){setcovariance <- 0}
  else {
    setcovariance <- 2*discr^2*acos(r/(2*discr)) - (r/2)*sqrt(4*discr^2-r^2)
  }
  return(setcovariance)
}
discr = 10
#plot to see what it looks like
plot(1:150/5,lapply(1:150/5,setcovdisc,discr=discr))
#is this right?
plot(1:150/5,2*discr^2*acos(-1:150/5/(2*discr)))
#this function looks right (at least if I've written it down correctly)
lines(1:150/5,((1:150/5)/2)*sqrt(4*discr^2-(1:150/5)^2))
#could plausibly cancel out the curviness of acos to make things straight close to 20
#compare to empirical overlaps of discs
areaintersecteddiscs <- function(displacement,discr){
  intersecteddiscs <- intersect.owin(disc(radius=10,centre=c(0,0)),disc(radius=10,centre=c(displacement,0)))
  return(area(intersecteddiscs))
}
plot(1:100/5,lapply(1:100/5,areaintersecteddiscs,discr=discr))
lines(1:150/5,lapply(1:150/5,setcovdisc,discr=discr))
#empirical matches theoretical formula! I'm happy with setcov discs functioning correctly :)

#combine into covariance of Boolean model
thcovDeterministicDiscs <- function(r,lambda,discr){
  expectedsetcovariance <- setcovdisc(r,discr)
  p <- 1-exp(-pi*discr^2*lambda)
  covariance <- 2*p-1+(1-p)^2*exp(lambda*expectedsetcovariance)
  return(covariance)
}
lambda <- 2.2064E-3
thcovDeterministicDiscs(10,lambda=lambda,discr=discr) #should be larger than 0.25 --> yes
thcovDeterministicDiscs(0,lambda=lambda,discr=discr) #which is p~0.5 :) when points get really close should tend towards the probability of the origin being in Xi
plot(1:110/5,lapply(1:110/5,thcovDeterministicDiscs,lambda = lambda, discr=discr))
#Matches: looks like the same curve as setcovdisc :), finishes at 0.25 ~ p^2, starts at p 

#theoretical covariance of a Boolean model is given (with proof) in Eqn 3.18 of Chiu 2014. Examples for particular grains can be found in B"ohm 2002
thcovDeterministicDiscs_vec <- function(X,Y,lambda,discr){
  rlist <- sqrt(X^2+Y^2)
  covar <- vector(length(rlist),mode="numeric")
  for (i in 1:length(rlist)){
    covar[i] <- thcovDeterministicDiscs(rlist[i],lambda=lambda,discr=discr)
  }
  return(covar)
}
#this plot should match the previous direction agnostic plot
lines(1:110/5,thcovDeterministicDiscs_vec(1:110/5,0,lambda=lambda,discr=discr))
#yes it does :)

#create a matrix map
xpts <- 0:250/10
ypts <- 0:250/10
mat <- outer(xpts,ypts,FUN="thcovDeterministicDiscs_vec",lambda=lambda,discr=discr) #rows correspond to xstep - just a quirk of outer!
image(mat) #rows in mat correspond to X, columns correspond to Y, but image() plots things s.t. rows increment along the X axis
#reflect out to all corners
mat <- mat[,c((ncol(mat)):2,1:ncol(mat))]
mat <- mat[c((nrow(mat)):2,1:nrow(mat)),]
image(mat)
#check that maximum is in the centre and that it is 0.5000069
max(mat)
arrayInd(which.max(mat),dim(mat))
#yes it is! :)

#compare to empirical covariance
#coverting theoretical to an image
thcov <- im(t(mat),xcol = c(-xpts[length(xpts):2],xpts),yrow=c(-ypts[length(ypts):2],ypts))
covar <- covarianceRACS(xi,w=Frame(xi))
plot(covar$covariance)
plot(covar$covariance,clipwin=Frame(thcov))
plot(thcov)
#looks mighty pixelated in comparison. I wonder how I could fix this?
#what with the stuff around 0.3 further away from the origin though?
plot(covar$covariance,clipwin=owin(xrange=c(-100,100),yrange=c(-100,100)))
#try to force higher resolution
covarHR <- covarianceRACS(as.mask(xi,dimyx=c(1000,1000)),w=Frame(xi))
plot(covarHR$covariance,clipwin=Frame(thcov))
plot(thcov)
plot(covarHR$covariance)
#it looks much nicer close to origin
#But the C(v)~0.3 stuff is in exactly the same spot as the low res version
#in fact can't tell the difference between plot(covarHR$covariance) and plot(covar$covariance)

#ok I'm happy that theorectical covariance is correct up to this point.
#convert into spectral domain to compare to spectral density
#note that in Bohm (2002) uses both a reduced covariance: 
#Cov(h) := P({o,h} in Xi) - p^2
#and the two-point Prob covariance C(h) := P({o,h} in Xi)
#the spectral density relationship applies to the former (Cov(h))
p <- 1-exp(-pi*discr^2*lambda)
thspecdens <- fft(as.matrix(thcov)- p^2)
image(Re(thspecdens))
image(Im(thspecdens))
image(abs(thspecdens)[1:15,1:15])
#according to Bohm 2002 just befor eqn 2.3 this should be a function that only maps to the Reals!
#but it is weird the imaginary part seems just as big - discussions with Myall - its just the phase of the signal - the peak isn't at 0 phase
#try padding with 0
M <- as.matrix(thcov)
nr <- nrow(M)
nc <- ncol(M)
thcovpad <- matrix(0, ncol=4*nc, nrow=4*nr)
thcovpad[1:nr, 1:nc] <- M
thspecdensPad <- fft(thcovpad)
filled.contour(Re(thspecdensPad))
image(Re(thspecdensPad)[1:20,1:20])
#there is slightly more structure now!
image(Im(thspecdensPad))
image(Im(thspecdensPad)[1:30,1:30])
image(abs(thspecdensPad)[1:30,1:30])
#looks very low res though. More padding == higher resolution cheat!
M <- as.matrix(thcov)
nr <- nrow(M)
nc <- ncol(M)
thcovpad <- matrix(0, ncol=8*nc, nrow=8*nr)
thcovpad[1:nr, 1:nc] <- M
thspecdensPad8 <- fft(thcovpad)
image(Re(thspecdensPad8)[1:20,1:20])
image(Im(thspecdensPad8)[1:20,1:20])
#but Imaginary element is still significant!
#Note this is an enormous matrix - 120MB!
#I think I need longer sample in spatial distance, means having a lower res sample
xptsLR <- 0:120/4
yptsLR <- 0:120/4
mat <- outer(xptsLR,yptsLR,FUN="thcovDeterministicDiscs_vec",lambda=lambda,discr=discr) #rows correspond to xstep - just a quirk of outer!
image(mat) #rows in mat correspond to X, columns correspond to Y, but image() plots things s.t. rows increment along the X axis
#reflect out to all corners
mat <- mat[,c((ncol(mat)):2,1:ncol(mat))]
mat <- mat[c((nrow(mat)):2,1:nrow(mat)),]
image(mat)
M <- mat-p^2
nr <- nrow(M)
nc <- ncol(M)
#pad with lots of 0!
thcovpad <- matrix(0, ncol=8*nc, nrow=8*nr)
thcovpad[1:nr, 1:nc] <- M
thspecdensPadLR <- fft(thcovpad)/(nr*8*nc*8)
image(Re(thspecdensPadLR))
image(Re(thspecdensPadLR)[1:30,1:30])
image(Im(thspecdensPadLR)[1:30,1:30])
#still imaginary part is about the same size. Both parts are oscillating too!!
#Their contours *should* be circular too :(

#for some reason the absolute value of theoretical spectral value is much curvier and nicer
image(abs(thspecdensPadLR)[1:30,1:30])
#xstep of thspecdensPadLR is 1/(xptsLR[1]-xptsLR[2])? But how does this resolution thing work?
plot(1:30,abs(thspecdensPadLR)[1,1:30])
#compare to smoothedspecdens
arrayInd(which.max(as.matrix(smoothedspecdens)),.dim = dim(smoothedspecdens))
plot(257:300,smoothedspecdens[258,257:300])
#see if can reproduce Bohm 2002 diagram
smoothedspecdens6 <-spectraldensity(xi,Frame(xi),6,dimyx=c(512,512))
smoothedspecdens8 <- spectraldensity(xi,Frame(xi),8,dimyx=c(512,512))
smoothedspecdens10 <- spectraldensity(xi,Frame(xi),10,dimyx=c(512,512))
lines(0:(300-257),smoothedspecdens6[258,257:300]/(6^2),col="red")
lines(0:(300-257),smoothedspecdens8[258,257:300]/(8^2),col="green")
plot(0:(300-257),smoothedspecdens10[258,257:300]/(10^2),type="l") #**note the scaling error in smoothed spectral density
max(smoothedspecdens10/(10^2))
max(abs(thspecdensPadLR))
#why the difference!?
lines(1:30,abs(thspecdensPadLR)[1,1:30]*max(smoothedspecdens10/(10^2))/max(abs(thspecdensPadLR)),col="blue")


dim(thspecdens)
plot(1:100,Re(thspecdens[1,1:100]))
plot(1:100,Im(thspecdens[1,1:100]))
nr = nrow(thspecdens)
nc = ncol(thspecdens)
thspecdens <- thspecdens[ ((-(nr-1)/2):((nr-1)/2)) %% (nr) + 1,((-nc/2):(nc/2)) %% (nc) + 1]
filled.contour(Re(thspecdens))