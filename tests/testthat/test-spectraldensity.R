context("Spectral Density")
#testspectraldensity

test_that("Spectal density calculations consistent with past calculations", {
#apply to heather
xi <- heather$coarse
specdens <- spectraldensity(xi,Frame(xi),20)
unsmspecdens <- unsmoothedspectraldensity(xi,Frame(xi))

#evaluate on a set of points
set.seed(31065461)
pttests <- rpoispp(0.025,win=owin(xrange=c(-10,10),yrange=c(-10,10)))

spdensamples <- specdens[pttests]
spdensamples_unsm <- unsmspecdens[pttests]
#saveRDS(spdensamples, "spectralden_sm_samples.RDS")
#saveRDS(spdensamples_unsm, "spectralden_unsm_samples.RDS")
expect_equal(spdensamples, readRDS("spectralden_sm_samples.RDS"))
expect_equal(spdensamples_unsm, readRDS("spectralden_unsm_samples.RDS"))
})
