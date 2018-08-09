context("Covariance Estimation")

  spatstat.options(npixel = 512) #to make default pixelisations higher resolution
  #estimate covariance of owin (covariance of image already estimated in help code)
  covarest.frowin.all <- balancedracscovariances(xi, obswin = w)

test_that("Covariance Estimation from an owin object matches the estimates from an im object", {
  
  covarest.frowin <- covarest.frowin.all[["pickaH"]]
  expect_is(covarest.frim, "im")
  expect_is(covarest.frowin, "im")

  covarest.diffs <- eval.im(a - b, harmonise.im(a = covarest.frowin, b = covarest.frim))

  expect_lt(mean(covarest.diffs), 0.1 * max(abs(covarest.frim)))
  expect_lt(max(abs(covarest.diffs)), 0.1 * max(abs(covarest.frim)))

})


#Test on a Boolean Model
#takes a few minutes
test_that("tradcovarest() matches theoretical covariance for Boolean Model", {

  #isotropise covar functions
  covarest.frim <- covarest.frim[owin(xrange = c(-1, 1) * 3.5 * discr, yrange = c(-1, 1) * 3.5 * discr), drop = TRUE]
  covarest.frowin.all <- lapply(covarest.frowin.all, 
				function(x) x[owin(xrange = c(-1, 1) * 3.5 * discr, yrange = c(-1, 1) * 3.5 * discr), drop = TRUE])
  covarest.frim.iso <- rotmean(covarest.frim, padzero = FALSE)
  covarest.frowin.all.iso <- lapply(covarest.frowin.all, rotmean, padzero = FALSE)

  truecovar.iso <- with.fv(covarest.frim.iso, bddcovar.iso(.x, lambda, discr),
                               fun = TRUE)
  names(truecovar.iso) <- c("r", "f")

  covarest.all.iso <- collapse.fv(c(list(truecovar = truecovar.iso),
				    list(frim = covarest.frim.iso),
				    covarest.frowin.all.iso), different = "f")

  #residual to theoretical covariance
  isocovarresid <- with.fv(covarest.all.iso, (. - truecovar) / truecovar, FUN = TRUE)
  #truncate residuals
  isocovarresid <- isocovarresid[isocovarresid$r < 3 * discr, ]

  #expect that this residual is smaller than 10% of the true covariance
  expect_lt(max(isocovarresid$pickaH),  0.1)
  expect_lt(max(isocovarresid$pickaintmult),  0.1)
  expect_lt(max(isocovarresid[, fvnames(isocovarresid, a = "."), drop = TRUE]), 0.1)
  #expect_lt(max(isocovarresid[, fvnames(isocovarresid, a = ".") != "none"]), 0.1)
  
})

  reset.spatstat.options()

test_that("unbalanced covariance estimation is symmetric for non-symmetric windows", {
  unsymmw <- union.owin(square(100),
	     shift.owin(square(100), vec = (c(100, 100))),
	     shift.owin(square(100), vec = (c(100, 200))),
	     shift.owin(square(100), vec = (c(200, 100)))
	     )
  
  setcovw <- setcov(unsymmw)
  expect_equal(max(abs(setcovw - reflect.im(setcovw))), 0)

  covarest.frowin <- balancedracscovariances(xi, obswin = unsymmw, modifications = "pickaH")[[1]]
  expect_equal(max(abs(covarest.frowin - reflect.im(covarest.frowin))), 0)
})

test_that("balancedracscovariances gives results in correct structure", {
  out <- balancedracscovariances(xi, obswin = w)
  expect_length(out, 8)
  expect_named(out)
  expect_s3_class(out, "imlist")
})

test_that("tradcovarest() errors properly", {
  lambda <- 4 * 2.2064E-3
  discr <- 5
  w <- owin(xrange = c(0, 100), yrange = c(0, 100))
  xi <- rbdd(lambda, discr, w)
  xiimg <- as.im(xi, W = w, eps = c(0.1, 0.1), na.replace = 0)
  xiimg[10, 10] <- 1.1
  expect_error(tradcovarest(xiimg), regexp = "Input xi has values other than 0, 1 or NA")
})
