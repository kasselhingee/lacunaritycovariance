# test twopt contag

library(stationaryracsinference, quietly = TRUE)

xi <- heather$coarse
covariance <- covariance(xi,Frame(xi))$covariance
twoptcontagion <- contagTwoPtProb(covariance)
p <- coveragefrac(xi,Frame(xi))
twoptcontagion <- contagTwoPtProb(covariance,p)
plot(twoptcontagion)
plot(twoptcontagion,clipwin=owin(xrange=c(-0.5,0.5),yrange=c(-0.5,0.5)),main="zoom")
v <- c(5,15)
twoptcontagion[ppp(v[1],v[2],window=owin(xrange=c(v[1]-1,v[1]+1),yrange=c(v[2]-1,v[2]+1)))]

contagTwoPtProb(covariance,p=p, v=c(5,15))
#result for both should be -1.384307