library(stationaryracsinference, quietly = TRUE)

#test units using a sine wave
A <- matrix(cos((1:1000)/10),nrow=10,ncol=1000,byrow=TRUE)
Aim <- im(A,xcol = (1:1000)/10, yrow = (1:10),
			unitname=c("metre","metres"))
fftA <- fft.im(Aim,padfactor=c(1,1))
#plot(Re(fftA),axes=TRUE,clipwin=owin(xrange=c(-4,4),yrange=c(-5,5)))
#The Fourier transform of cos(x) is infinite at theta = 1 (because integrating cos^2(x) out to infinity),
# and 0 elsewhere. So result should be an image with units 1/m, with a single large peak at x=1 
maxindex <- which(as.matrix(Re(fftA)==max(Re(fftA))), arr.ind=TRUE)
abs(fftA$xcol[maxindex[2]]+1)<fftA$xstep
abs(fftA$yrow[maxindex[1]]-0)<fftA$ystep


#test height using a Gaussian
gauss <- function(x,coeff=1){ return(exp(-coeff*x*x)) }

B <- matrix(gauss((-1000:1000)/10,coeff=4),nrow=10,ncol=2001,byrow=TRUE)
Bim <- im(B,xcol = (-1000:1000)/10, yrow = 1:10,
			unitname=c("metre","metres"))
fftB <- fft.im(Bim,padfactor=c(1,1))
#the result should approximate sqrt(pi/coeff)*exp(-t^2/4coeff) along the x-axis
#At t=0 it is sqrt(pi/4)=0.8862269. However because the Gaussian is repeated 10 times (10 rows) 
#the height should be 10*0.8862269
abs(max(Re(fftB)) - 10*sqrt(pi/4)) < 1E-8