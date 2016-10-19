#testspectraldensity

suppressPackageStartupMessages(library(stationaryracsinference))

#apply to heather
xi <- heather$coarse
specdens <- spectraldensity(xi,Frame(xi),20)
unsmspecdens <- unsmoothedspectraldensity(xi,Frame(xi))

#evaluate on a set of points
set.seed(31065461)
pttests <- rpoispp(0.025,win=owin(xrange=c(-10,10),yrange=c(-10,10)))

specdens[pttests]
unsmspecdens[pttests]
