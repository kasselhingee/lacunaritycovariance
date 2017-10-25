
lambda <- 4 * 2.2064E-3
discr <- 5
w <- owin(xrange = c(0, 100) * 3, yrange = c(0, 100) * 3)
xi <- rbdd(lambda, discr, w)
xiimg <- as.im(xi, W = w, eps = c(0.1, 0.1), na.replace = 0)
xi.p <- sum(xiimg) / sum(is.finite(xiimg$v))
#estimate covariance
spatstat.options(npixel = 512) #to make default pixelisations higher resolution
covarest.frim <- racscovariance(xiimg)
reset.spatstat.options()