
#Test on a Boolean Model
#takes a few minutes
test_that("covariance() matches theoretical covariance for Boolean Model", {
  lambda <- 4 * 2.2064E-3
  discr <- 5
  w <- owin(xrange = c(0, 100) * 2, yrange = c(0, 100) * 2)
  xi <- rbdd(lambda, discr,w)
  xiimg <- as.im(xi, W = w, eps = c(0.1, 0.1), na.replace = 0)
  #estimate covariance
  spatstat.options(npixel = 512) #to make default pixelisations higher resolution
  covarest.frowin <- covariance(xi, obswin = w)
  covarest.frim <- covariance(xiimg)
  
  expect_is(covarest.frim, "im")
  expect_is(covarest.frowin, "im")
  
  covarest.diffs <- covarest.frowin - covarest.frim
  
  expect_lt(mean(covarest.diffs), 0.1 * max(abs(covarest.frim)))
  expect_lt(max(abs(covarest.diffs)), 0.1 * max(abs(covarest.frim)))
  
  #isotropise above functions
  covarest.frim.iso <- rotmean(covarest.frim)
  covarest.frowin.iso <- rotmean(covarest.frowin)

  truecovar.iso <- with.fv(covarest.frim.iso, bdd_covar.iso(.x, lambda, discr),
                               fun = TRUE)
  
  #residual to theoretical covariance
  isocovarresid.frim <- eval.fv((covarest.frim.iso - truecovar.iso) / truecovar.iso, equiv = c(f = "col1"))
  isocovarresid.frowin <- eval.fv((covarest.frowin.iso - truecovar.iso) / truecovar.iso, equiv = c(f = "col1"))
  
  
  #expect that this residual is smaller than 10% of the true covariance
  maxisocovarresid.frim <- max(eval.fv(abs(isocovarresid.frim))[isocovarresid.frim$r < 3 * discr,])
  expect_lt(maxisocovarresid.frim,  0.1)
  
  maxisocovarresid.frowin <- max(eval.fv(abs(isocovarresid.frowin))[isocovarresid.frowin$r < 3 * discr,])
  expect_lt(maxisocovarresid.frowin, 0.1)
  
  reset.spatstat.options()
})
