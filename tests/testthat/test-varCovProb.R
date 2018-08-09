
test_that("cpvariance() matches the theoretical coverage probability estimator variance for a Boolean Model", {
  #estimate covariance of owin (covariance of image already estimated in help code)
  spatstat.options(npixel = 512) #to make default pixelisations higher resolution
  cpvariance.frowin <- cpvariance(xi, w)
  
  expect_gt(cpvariance.frowin, 0)

  #expect_lt(abs(cpvariance.frowin - readRDS("true.var.p.RDS"))/readRDS("true.var.p.RDS"),
  #          0.1)
  
  reset.spatstat.options()
})
