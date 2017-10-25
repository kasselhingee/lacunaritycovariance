
#Test on a Boolean Model
#takes a few minutes
test_that("racscovariance() matches theoretical covariance for Boolean Model", {
  #estimate covariance of owin (covariance of image already estimated in help code)
  spatstat.options(npixel = 512) #to make default pixelisations higher resolution
  covarest.frowin <- racscovariance(xi, obswin = w)

  expect_is(covarest.frim, "im")
  expect_is(covarest.frowin, "im")

  covarest.diffs <- eval.im(a - b, harmonise.im(a = covarest.frowin, b = covarest.frim))

  expect_lt(mean(covarest.diffs), 0.1 * max(abs(covarest.frim)))
  expect_lt(max(abs(covarest.diffs)), 0.1 * max(abs(covarest.frim)))

  #isotropise above functions
  covarest.frim <- covarest.frim[owin(xrange = c(-1, 1) * 3.5 * discr, yrange = c(-1, 1) * 3.5 * discr), drop = TRUE]
  covarest.frowin <- covarest.frowin[owin(xrange = c(-1, 1) * 3.5 * discr, yrange = c(-1, 1) * 3.5 * discr), drop = TRUE]
  covarest.frim.iso <- rotmean(covarest.frim, padzero = FALSE)
  covarest.frowin.iso <- rotmean(covarest.frowin, padzero = FALSE)

  truecovar.iso <- with.fv(covarest.frim.iso, bddcovar.iso(.x, lambda, discr),
                               fun = TRUE)

  #residual to theoretical covariance
  isocovarresid.frim <- eval.fv( (covarest.frim.iso - truecovar.iso) / truecovar.iso, equiv = c(f = "col1"))
  isocovarresid.frowin <- eval.fv( (a - b) / b, envir = harmonise.fv(a = covarest.frowin.iso, b = truecovar.iso), equiv = c(f = "col1"))


  #expect that this residual is smaller than 10% of the true covariance
  maxisocovarresid.frim <- max(eval.fv(abs(isocovarresid.frim))[isocovarresid.frim$r < 3 * discr, ])
  expect_lt(maxisocovarresid.frim,  0.1)
  
  maxisocovarresid.frowin <- max(eval.fv(abs(isocovarresid.frowin))[isocovarresid.frowin$r < 3 * discr, ])
  expect_lt(maxisocovarresid.frowin, 0.1)

  reset.spatstat.options()
})

test_that("racscovariance() errors properly", {
  lambda <- 4 * 2.2064E-3
  discr <- 5
  w <- owin(xrange = c(0, 100), yrange = c(0, 100))
  xi <- rbdd(lambda, discr, w)
  xiimg <- as.im(xi, W = w, eps = c(0.1, 0.1), na.replace = 0)
  xiimg[10, 10] <- 1.1
  expect_error(racscovariance(xiimg), regexp = "Input xi has values other than 0, 1 or NA")
})
