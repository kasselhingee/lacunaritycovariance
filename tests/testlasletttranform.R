# laslett transform script

library(stationaryracsinference, quietly = TRUE)

#test on a known true boolean model
set.seed(135498156)
xi <- rBooleanDetermDiscs(2.2064E-3,10,owin(xrange=c(0,500),yrange=c(0,200)))
xim <- as.mask(xi,eps=c(1,1))
lltp <- findlowlefttangentpts(xi)
length(lltp$X)
lltp
ltxi <- lasletttransform(xi)
ltxi$n
Ks <- Kest(ltxi)
max(abs(Ks$theo-Ks$trans))
