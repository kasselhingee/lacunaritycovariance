context("Estimation - second order properties")

spatstat.options(npixel = 2^3)
xi <- heather$coarse
xiim_verytoy <- as.im(xi, value = TRUE, na.replace = FALSE, eps = 2)

test_that("secondorderprop() produces output of correct class with all estimators", {
  skip_on_cran()  #takes 12 seconds on Kassel's laptop
  expect_warning(secondests <- secondorderprops(xiim_verytoy, 
    gblargs = list(boxwidths = seq(1, 10, by = 0.1)),
    covarargs = list(estimators = "all"),
    cencovarargs = list(estimators = "all"),
    paircorrargs = list(estimators = "all"),
    returnrotmean = FALSE),
    regexp = "harmonised")
  expect_s3_class(secondests$GBL, "fv")
  expect_s3_class(secondests$covariance, "imlist")
  expect_s3_class(secondests$cencovariance, "imlist")
  expect_s3_class(secondests$paircorr, "imlist")
  expect_gt(length(secondests$GBL), 3)
  expect_gt(length(secondests$covariance), 3)
  expect_gt(length(secondests$cencovariance), 3)
  expect_gt(length(secondests$paircorr), 3)
})

test_that("secondorderprop() produces empty list with full defaults", { 
  secondests <- secondorderprops(xiim_verytoy)
  expect_is(secondests, "list")
  expect_length(secondests, 0)
})

test_that("secondorderprop() produces of fv objects for rotmean with multiple estimators", {
  skip_on_cran() #takes 13 seconds on Kassel's laptop
  expect_warning(secondests <- secondorderprops(xiim_verytoy, 
    gblargs = list(boxwidths = seq(1, 10, by = 0.1)),
    covarargs = list(estimators = "all"),
    cencovarargs = list(estimators = "all"),
    paircorrargs = list(estimators = "all"),
    returnrotmean = TRUE),
    regexp = "harmonised")
  expect_s3_class(secondests$GBL, "fv")
  expect_s3_class(secondests$covariance, "fv")
  expect_s3_class(secondests$cencovariance, "fv")
  expect_s3_class(secondests$paircorr, "fv")
  expect_gt(length(secondests$GBL), 3)
  expect_gt(length(secondests$covariance), 3)
  expect_gt(length(secondests$cencovariance), 3)
  expect_gt(length(secondests$paircorr), 3)
})
  
test_that("secondorderprop() produces of fv objects for rotmean with single estimators", { 
  #test with rotmean and single estimators
  secondests <- secondorderprops(xiim_verytoy, 
    gblargs = list(boxwidths = seq(1, 10, by = 0.1), estimators = "GBLcc.pickaH"),
    covarargs = list(estimators = "pickaH"),
    cencovarargs = list(estimators = "pickaH"),
    paircorrargs = list(estimators = "pickaH"),
    returnrotmean = TRUE)
  expect_s3_class(secondests$GBL, "fv")
  expect_s3_class(secondests$covariance, "fv")
  expect_s3_class(secondests$cencovariance, "fv")
  expect_s3_class(secondests$paircorr, "fv")
  expect_length(secondests$GBL, 2)
  expect_length(secondests$covariance, 2)
  expect_length(secondests$cencovariance, 2)
  expect_length(secondests$paircorr, 2)
})
  
test_that("secondorderprop() single images without rotmean and when only one estimator requested", { 
  skip_on_cran()  #non-important function doesn't need testing on CRAN
  secondests <- secondorderprops(xiim_verytoy, 
    gblargs = list(boxwidths = seq(1, 10, by = 0.1), estimators = "GBLcc.pickaH"),
    covarargs = list(estimators = "pickaH"),
    cencovarargs = list(estimators = "pickaH"),
    paircorrargs = list(estimators = "pickaH"),
    returnrotmean = FALSE)
  expect_s3_class(secondests$GBL, "fv")
  expect_s3_class(secondests$covariance, "imlist")
  expect_s3_class(secondests$cencovariance, "imlist")
  expect_s3_class(secondests$paircorr, "imlist")
  expect_length(secondests$GBL, 2)
  expect_length(secondests$covariance, 1)
  expect_length(secondests$cencovariance, 1)
  expect_length(secondests$paircorr, 1)
})

test_that("secondorderprop() gives single fv object when only GBLemp", {
  skip_on_cran()  #non-important function doesn't need testing on CRAN
  secondests <- secondorderprops(xiim_verytoy, 
    gblargs = list(boxwidths = seq(1, 10, by = 0.1), estimators = "GBLemp"))
  expect_s3_class(secondests$GBL, "fv")
  expect_length(secondests, 1)
})

test_that("secondorderprop() gives data frame for discs", {
  skip_on_cran()  #non-important function doesn't need testing on CRAN
  discs <- lapply(seq(1, 10, by = 0.5), function(x) disc(r = x))
  secondests <- secondorderprops(xiim_verytoy, 
    gblargs = list(boxwidths = discs, estimators = "GBLc"))
  expect_s3_class(secondests$GBL, "data.frame")
  expect_length(secondests, 1)
})

  test_that("secondorderprop() gives data frame for discs, multiple estimators", {
  skip_on_cran()  #non-important function doesn't need testing on CRAN
  discs <- lapply(seq(1, 10, by = 0.5), function(x) disc(r = x))
  secondests <- secondorderprops(xiim_verytoy, 
    gblargs = list(boxwidths = discs, estimators = c("GBLc", "GBLcc.pickaH", "GBLcc.mattfeldt", "GBLcc.pickaint", "GBLg.pickaint", "GBLg.mattfeldt")))
  expect_s3_class(secondests$GBL, "data.frame")
  expect_equal(nrow(secondests$GBL), length(discs))
})
  
  test_that("secondorderprop() gives error when discs and GBLemp", {
  discs <- lapply(seq(1, 10, by = 0.5), function(x) disc(r = x))
  expect_error(secondorderprops(xiim_verytoy, 
    gblargs = list(boxwidths = discs, estimators = "GBLemp")),
    regexp = "non-numeric")
})

  test_that("secondorderprop() ignores centred covariance when requested to", {
  skip_on_cran()  #non-important function doesn't need testing on CRAN
  secondests <- secondorderprops(xiim_verytoy, 
    gblargs = list(boxwidths = seq(1, 10, by = 0.1), estimators = "GBLcc.pickaH"),
    covarargs = list(estimators = "pickaH"),
    paircorrargs = list(estimators = "pickaH"),
    returnrotmean = FALSE)
  expect_s3_class(secondests$GBL, "fv")
  expect_s3_class(secondests$covariance, "imlist")
  expect_null(secondests$cencovariance)
  expect_s3_class(secondests$paircorr, "imlist")
  expect_length(secondests$GBL, 2)
  expect_length(secondests$covariance, 1)
  expect_length(secondests$paircorr, 1)
})

reset.spatstat.options()
