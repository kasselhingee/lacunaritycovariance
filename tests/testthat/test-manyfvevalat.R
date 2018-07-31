context("Many fv eval")

test_that("manyfvevalat gives the values of the indivual functions", {
fv1 <- Hest(heather$coarse)
fv2 <- Hest(complement.owin(heather$coarse))
fvlist <- list(fv1,fv2) 

fvfunc01 <- as.function.fv(fv1, value="km")
fvfunc02 <- as.function.fv(fv2, value="km")

expect_equal(manyfvevalat(fvlist,0.19,value="km"),
	     list(fvfunc01(0.19), fvfunc02(0.19)))
})
