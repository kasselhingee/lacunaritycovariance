
lambda <- 4 * 2.2064E-3
discr <- 5
w <- owin(xrange = c(0, 100), c(0, 100))
if (identical(Sys.getenv("NOT_CRAN"), "true")){
  w <- owin(xrange = c(0, 100) * 3, yrange = c(0, 100) * 3)
}
xi <- rbdd(lambda, discr, w)
if (identical(Sys.getenv("NOT_CRAN"), "true")){
    xiimg <- as.im(xi, W = w, eps = c(0.5, 0.5), na.replace = 0)
} else {
    xiimg <- as.im(xi, W = w, eps = c(2, 2), na.replace = 0)
}
xi.p <- sum(xiimg) / sum(is.finite(xiimg$v))
#estimate covariance
covarest.frim <- racscovariance(xiimg, estimators = "pickaH", drop = TRUE)

# # # saved calculation of true coverage fraction variance
# true.covar <- bddcovar(covarest.frowin$xrange,
#         covarest.frowin$yrange,
#         eps = c(covarest.frowin$xstep/2, covarest.frowin$ystep/2),
#         lambda = lambda,
#         discr = discr)
# true.p <- bddcoverageprob(lambda, discr)
# 
# setcovB <- setcov(w)
# integrand <- eval.im((thcovar-(true.p^2))*setcovB, harmonize = TRUE)
# true.var.p <- ((1/(area.owin(w))^2)*sum(integrand)*integrand$xstep*integrand$ystep)
# saveRDS(true.var.p, "true.var.p.RDS")
