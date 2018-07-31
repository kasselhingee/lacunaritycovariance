context("Contagion Computations")
# test twopt contag

test_that("Two Point (Covariance) Contagion is consitent with past computations", {
xi <- heather$coarse
covariance <- racscovariance(xi,Frame(xi))
twoptcontagion <- contagTwoPtProb(covariance)
p <- coveragefrac(xi,Frame(xi))
twoptcontagion <- contagTwoPtProb(covariance,p)
#plot(twoptcontagion)
#plot(twoptcontagion,clipwin=owin(xrange=c(-0.5,0.5),yrange=c(-0.5,0.5)),main="zoom")
v <- c(5,15)
atv_1 <- twoptcontagion[ppp(v[1],v[2],window=owin(xrange=c(v[1]-1,v[1]+1),yrange=c(v[2]-1,v[2]+1)))]

atv_2 <- contagTwoPtProb(covariance,p=p, v=c(5,15))
#result for both should be -1.384307
expect_equal(atv_1, -1.384307, tolerance = 1E-3)
expect_equal(atv_2, -1.384307, tolerance = 1E-3)
})
