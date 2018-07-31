context("Laslett Transform")
# laslett transform script

test_that("Laslett Transform Matches Past Results for Heather Data", {
#test on a known true boolean model
xi <- heather$fine
xim <- as.mask(xi)
lltp <- findlowlefttangentpts(xi)
expect_length(lltp$X, 95)
})

test_that("Laslett Transform produces a Poisson process from the Boolean model", {
#check that Boolean model gives a Poisson point process according to K function estimates
xi <- rbdd(2.2064E-3,10,owin(xrange=c(0,500),yrange=c(0,200)))
xim <- as.mask(xi,eps=c(1,1))
ltxi <- lasletttransform(xi)
#if I was doing this well I would use a hypothesis test from spatstat:
mdtestout <- mad.test(ltxi, Lest, verbose = FALSE)
expect_gt(mdtestout$p.value, 0.01) #this test should only fail 1 in every 100
})
