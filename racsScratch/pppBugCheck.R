#checking a ppp creation bug found when writing Laslett Transform
#actually it was me needing to put window = w in the argument to the function

library(devtools)
install_github('spatstat/spatstat')
library(spatstat)

#create a list of points
pttable <- data.frame(X = 1:5,Y = 1:5)

w <- owin(xrange = c(0,6),yrange=c(0,6))

pointPat <- ppp(pttable$X,pttable$Y,w)
#Error in owin(...) : 
#  If one of xrange, yrange is specified then both must be.

pointPat <- ppp(pttable$X,pttable$Y,window=w)
#this works without problem!
