suppressPackageStartupMessages(library(stationaryracsinference))

covar <- covariance(heather$coarse)
p <- area(heather$coarse)/area(Frame(heather$coarse))
sidelengths <- 2.2
lac <- lac(sidelengths,covar,p)
lac$MVL

#compare to calculation from spatstat's setcov:
abs(lac$MVL - mvlc(list(square(2.2)),covar,p) ) < 0.01

#' #Test on a Boolean Model #takes a few minutes
lambda <- 4*2.2064E-3
discr <- 5
w <- owin(xrange=c(0,100)*2,yrange=c(0,100)*2)
xi <- rBooleanDetermDiscs(lambda,discr,w)
xiimg <- as.im(xi, W=w, eps=c(0.1,0.1), na.replace=0)
#estimate lacunarity
sidelengths <- seq(0.5,50,by=1)
mvl.est <- mvlc(sidelengths,xiim=xiimg)
#theoretical lacunarity very different because window is small **I think
thcovariance <- thcovarDeterministicDiscs(
                 xrange=c(-10,10)*4,
	    yrange=c(-10,10)*4,
	    eps=c(0.1,0.1),lambda,discr)
#the following makes sure there is a border of NA values so that mvlc knows when it is using something outside the image
thcovariance[setminus.owin(Frame(thcovariance),erosion(Frame(thcovariance),0.1))] <- NA
thcoverageprob <- booldetermdiscs_truecoveragefrac(lambda,discr)
mvl.th <- mvlc(sidelengths,thcovariance,thcoverageprob)
abs(max(eval.fv(mvl.est-mvl.th)))  < 0.05
