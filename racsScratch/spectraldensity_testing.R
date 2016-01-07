#testing spectral density calculations
library(stationaryracsinference)

#following on partially through checking_spectraldensity_calcs.Rmd
#after generating a low resolution version of the theoretical covariance
#mat is this low res version

image(mat)
M <- mat-p^2
nr <- nrow(M)
nc <- ncol(M)
#pad with lots of 0!
thcovpad <- matrix(0, ncol=8*nc, nrow=8*nr)
thcovpad[1:nr, 1:nc] <- M
#calculate the scaling factor (=distance between (grid points)^2) - from my note book p95
scalefactorX <- (xptsLR[2]-xptsLR[1])
scalefactorY <- (yptsLR[2]-yptsLR[1]) 
thspecdensPadLR <- scalefactorX*scalefactorY*fft(thcovpad)
plot(1:30,abs(thspecdensPadLR)[1,1:30]) #this is the same scale as Bohm2002 plots!


#testing phase shift ... it didn't work
spectralpointsX <- (2*pi/(scalefactorX*nrow(thcovpad))) * (0:(nrow(thcovpad)-1))
spectralpointsY <- (2*pi/(scalefactorY*ncol(thcovpad))) * (0:(ncol(thcovpad)-1))
firstpointX <- -xptsLR[length(xptsLR)]
firstpointY <- -yptsLR[length(yptsLR)]

phaseshiftmatrix <- exp(1i*outer(spectralpointsX*firstpointX,spectralpointsY*firstpointY,FUN = "+"))
filled.contour(Re(shiftmatrix[1:30,1:30]))
filled.contour(Im(shiftmatrix[1:30,1:30]))

filled.contour(Re(thspecdensPadLR[1:30,1:30]))
filled.contour(Im(thspecdensPadLR[1:30,1:30]))
image(Im(thspecdensPadLR[1:30,1:30]))

#but appling phase shiftvector doesn't shift it into the reals :(
image(Re(thspecdensPadLR[1:30,1:30]*shiftmatrix[1:30,1:30]))
image(Im(thspecdensPadLR[1:30,1:30]*shiftmatrix[1:30,1:30]))

plot(thspecdensPadLR[1,30])
lines(shiftmatrix[1,30])

unsmsd <- unsmoothedspectraldensity(xi,Frame(xi),dimyx=c(512,512))
plot(unsmsd,clipwin=owin(xrange = c(-0.1,0.1),yrange=c(-0.1,0.1)))
plot(0:(300-257),unsmsd[258,257:300]) #this looks like the right sort of scale

#test smoothed version. But it had memory issues - the bandwidth is waaay too large
smoothedsd01 <-spectraldensity(xi,Frame(xi),0.01,dimyx=c(512,512))
smoothedsd02 <-spectraldensity(xi,Frame(xi),0.02,dimyx=c(512,512))
smoothedsd03 <-spectraldensity(xi,Frame(xi),0.03,dimyx=c(512,512))
smoothedsd05 <-spectraldensity(xi,Frame(xi),0.05,dimyx=c(512,512))
smoothedsd10 <-spectraldensity(xi,Frame(xi),0.10,dimyx=c(512,512)) #this bandwidth is massive - much larger than the spike
smoothedsd20 <-spectraldensity(xi,Frame(xi),0.20,dimyx=c(512,512))
plot(smoothedsd01,clipwin=owin(xrange = c(-0.1,0.1),yrange=c(-0.1,0.1)))
plot(smoothedsd10,clipwin=owin(xrange = c(-0.1,0.1),yrange=c(-0.1,0.1)))

#compare sizes

plot(0:(300-257),smoothedsd01[258,257:300],type="l",col="red") #**note the scaling error in smoothed spectral density
lines(0:(300-257),smoothedsd02[258,257:300],col="pink") #**note the scaling error in smoothed spectral density
lines(0:(300-257),smoothedsd03[258,257:300],col="yellow") #**note the scaling error in smoothed spectral density
lines(0:(300-257),smoothedsd05[258,257:300],col="green") #**note the scaling error in smoothed spectral density
lines(0:(300-257),smoothedsd10[258,257:300],col="blue") #**note the scaling error in smoothed spectral density

lines(1:30,abs(thspecdensPadLR)[1,1:30]*max(smoothedsd03)/max(abs(thspecdensPadLR)),col="black",lwd=2)

#something is weird with scaling of smoother
#Test on unit image
flatSignal=as.im(owin(xrange = c(-1,1),yrange=c(-1,1)),dimyx=c(100,100))
bandwidth = 0.1
xstep = flatSignal$xstep
ystep = flatSignal$ystep
supportwidth = 3*bandwidth
X <- seq(0,supportwidth*1.5+xstep,by=xstep) #much larger than support width to avoid boundary issues?
Y <- seq(0,supportwidth*1.5+ystep,by=ystep)
  mat <- outer(X/bandwidth,Y/bandwidth,FUN="EpanechnikovFcn") #rows correspond to xstep - just a quirk of outer!
  mat <- mat/(bandwidth^2) #to account for the scaling of the kernel - so that it all adds to 1
  mat <- t(mat) #columns correspond to changes in X, rows correspond to changes in Y!
  #reflect out to all corners
  mat <- mat[,c((ncol(mat)):2,1:ncol(mat))]
  mat <- mat[c((nrow(mat)):2,1:nrow(mat)),]
  kernelfcn <- im(mat,xcol=c(-X[length(X):2],X),yrow=c(-Y[length(Y):2],Y))
#apply kernel using convolve.im
smoothedsig <- convolve.im(flatSignal,kernelfcn)
plot(smoothedsig)
max(smoothedsig)
smoothedsig <- smspecdens[Frame(specdens)]

#seems it was my dodgy integration and weight kernel function incorrectly
detach(package:stationaryracsinference,unload=TRUE)
library(stationaryracsinference)
#test smoothed version. But it had memory issues - the bandwidth is waaay too large
smoothedsd01 <-spectraldensity(xi,Frame(xi),0.01,dimyx=c(512,512))
smoothedsd015 <-spectraldensity(xi,Frame(xi),0.015,dimyx=c(512,512))
smoothedsd02 <-spectraldensity(xi,Frame(xi),0.02,dimyx=c(512,512))
smoothedsd03 <-spectraldensity(xi,Frame(xi),0.03,dimyx=c(512,512))
smoothedsd05 <-spectraldensity(xi,Frame(xi),0.05,dimyx=c(512,512))
smoothedsd10 <-spectraldensity(xi,Frame(xi),0.10,dimyx=c(512,512)) #this bandwidth is massive - much larger than the spike
smoothedsd20 <-spectraldensity(xi,Frame(xi),0.20,dimyx=c(512,512))
plot(smoothedsd01,clipwin=owin(xrange = c(-0.1,0.1),yrange=c(-0.1,0.1)))
plot(smoothedsd10,clipwin=owin(xrange = c(-0.1,0.1),yrange=c(-0.1,0.1)))

#compare sizes
plot(0:(300-257),unsmsd[258,257:300])
plot(1:30,abs(thspecdensPadLR)[1,1:30],col="black",lwd=2)
lines(0:(300-257),smoothedsd01[258,257:300],type="l",col="red") #**note the scaling error in smoothed spectral density
lines(0:(300-257),smoothedsd015[258,257:300],col="pink") #**note the scaling error in smoothed spectral density
lines(0:(300-257),smoothedsd03[258,257:300],col="yellow") #**note the scaling error in smoothed spectral density
lines(0:(300-257),smoothedsd05[258,257:300],col="green") #**note the scaling error in smoothed spectral density
lines(0:(300-257),smoothedsd10[258,257:300],col="blue") #**note the scaling error in smoothed spectral density

spectraldensity() still doesn't seem to be a good estimate of spectral density!!