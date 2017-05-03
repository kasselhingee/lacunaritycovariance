# test integratebyradius
suppressPackageStartupMessages(library(stationaryracsinference))

testim <- square(r=5)
testim <- setminus.owin(testim,square(r=2))
#plot(testim,axes=TRUE)
kfcn <- integratebyradius(c(2,2),as.im(testim,eps=0.1,na.replace=0))
#plot(kfcn)
#lines(kfcn$r,0.75*pi*kfcn$r*kfcn$r,col="red",lty="dashed")
max(abs(kfcn$y - 0.75*pi*kfcn$r*kfcn$r))

kfcn <- integratebyradius(c(2,2),
                  as.im(testim,eps=c(0.1,0.5),na.replace=0))
#plot(kfcn)
#lines(kfcn$r,0.75*pi*kfcn$r*kfcn$r,col="red",lty="dashed")
max(abs(kfcn$y - 0.75*pi*kfcn$r*kfcn$r))