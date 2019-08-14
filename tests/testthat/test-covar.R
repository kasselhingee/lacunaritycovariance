context("Estimation - Covariance")

  if (identical(Sys.getenv("NOT_CRAN"), "true")){
    spatstat.options(npixel = 512) #to make default pixelisations higher resolution
    proptol = 0.1
  } else {
    spatstat.options(npixel = xiimg$dim[[1]])  #matches the number of pixels of the discretised xi
    proptol = 0.1
  }
  #estimate covariance of owin (covariance of image already estimated in help code)
  covarest.frowin.all <- racscovariance(xi, obswin = w)

test_that("Covariance Estimation from an owin object matches the estimates from an im object", {
  
  covarest.frowin <- covarest.frowin.all[["pickaH"]]
  expect_is(covarest.frim, "im")
  expect_is(covarest.frowin, "im")

  covarest.diffs <- eval.im(a - b, harmonise.im(a = covarest.frowin, b = covarest.frim))

  expect_lt(mean(covarest.diffs) / max(abs(covarest.frim)), proptol)
  expect_lt(max(abs(covarest.diffs)) / max(abs(covarest.frim)), proptol) 

})


#Test on a Boolean Model
#takes a few minutes
test_that("racscovariance pickaH method matches theoretical covariance for Boolean Model", {
  skip_on_cran()

  #isotropise covar functions
  covarest.frim <- covarest.frim[owin(xrange = c(-1, 1) * 3.5 * discr, yrange = c(-1, 1) * 3.5 * discr), drop = TRUE]
  covarest.frowin.all <- lapply(covarest.frowin.all, 
				function(x) x[owin(xrange = c(-1, 1) * 3.5 * discr, yrange = c(-1, 1) * 3.5 * discr), drop = TRUE])
  covarest.frim.iso <- rotmean(covarest.frim, padzero = FALSE)
  covarest.frowin.all.iso <- lapply(covarest.frowin.all, rotmean, padzero = FALSE)

  truecovar.iso <- with.fv(covarest.frim.iso, bddcovar.iso(.x, lambda, discr),
                               fun = TRUE)
  names(truecovar.iso) <- c("r", "f")

  covarest.all.iso <- collapse.fv(harmonise.fv(c(list(truecovar = truecovar.iso),
				    list(frim = covarest.frim.iso),
				    covarest.frowin.all.iso)), different = "f")

  #residual to theoretical covariance
  isocovarresid <- with.fv(covarest.all.iso, (. - truecovar) / truecovar, FUN = TRUE)
  #truncate residuals
  isocovarresid <- isocovarresid[isocovarresid$r < max(isocovarresid$r/2), ]  #the condition removes large radii relative to window size

  #expect that this residual is smaller than 10% of the true covariance
  expect_lt(max(isocovarresid$pickaH),  proptol)
  expect_lt(max(isocovarresid$pickaint),  proptol)
  expect_lt(max(isocovarresid$mattfeldt),  proptol)
  expect_lt(max(isocovarresid[, fvnames(isocovarresid, a = "."), drop = TRUE]), proptol)
  #expect_lt(max(isocovarresid[, fvnames(isocovarresid, a = ".") != "plugin"]), 0.1)
  
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

  covarest.frowin <- racscovariance(xi, obswin = unsymmw, estimators = "pickaH", drop = TRUE)
  expect_equal(max(abs(covarest.frowin - reflect.im(covarest.frowin))), 0)
})

test_that("racscovariance gives results in correct structure", {
  out <- racscovariance(xi, obswin = w)
  expect_length(out, 4)
  expect_named(out)
  expect_s3_class(out, "imlist")
})

test_that("racscovariance errors correctly", {
  expect_error(racscovariance(heather$coarse), regexp = "owin")
  expect_error(racscovariance(as.im(heather$coarse, value = TRUE, na.replace = FALSE),
                 obswin = shift.owin(disc(), vec = c(3, 3))))
})

test_that("plugincvc errors when not passed a binary map", {
  xi <- heather$coarse
  xiimg <- as.im(xi, W = w, eps = c(0.1, 0.1), na.replace = 0)
  xiimg[5, 10] <- 1.1
  expect_error(plugincvc(xiimg), regexp = "Input xi has values other than 0, 1 or NA")
})

test_that("plugincvc errors when passed a foreground set and no observation window", {
  xi <- heather$coarse
  expect_error(plugincvc(xi), regexp = "obswin")
})

test_that("plugincvc maintains unitnames", {
  xi <- heather$coarse
  cvcest <- plugincvc(xi, obswin = Frame(xi))
  expect_equal(unitname(xi), unitname(cvcest))
})
