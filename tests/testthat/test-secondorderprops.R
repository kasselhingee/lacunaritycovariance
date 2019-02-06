context("Estimation - second order properties")

test_that("secondorderprop() produces output of the expected classes given different arguments", {
  xi <- heather$coarse
  xiim <- as.im(xi, value = TRUE, na.replace = FALSE)

  #operates with all props going "everythingmode"
  expect_warning(secondests <- secondorderprops(xiim, 
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

  #test with full defaults
  secondests <- secondorderprops(xiim)
  expect_is(secondests, "list")
  expect_length(secondests, 0)

  #test with rotmean and multiple estimators
  expect_warning(secondests <- secondorderprops(xiim, 
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
  
  #test with rotmean and single estimators
  secondests <- secondorderprops(xiim, 
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
  
  #test without rotmean and single estimators
  secondests <- secondorderprops(xiim, 
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

  #test with only GBLemp and partially NULL for the other ones
  secondests <- secondorderprops(xiim, 
    gblargs = list(boxwidths = seq(1, 10, by = 0.1), estimators = "GBLemp"))
  expect_s3_class(secondests$GBL, "fv")
  expect_length(secondests, 1)
  
  #test with discs and no GBLemp
  discs <- lapply(seq(1, 10, by = 0.5), function(x) disc(r = x))
  secondests <- secondorderprops(xiim, 
    gblargs = list(boxwidths = discs, estimators = "GBLc"))
  expect_s3_class(secondests$GBL, "data.frame")
  expect_length(secondests, 1)
  
  #test with discs and multpile estimators,and no GBLemp
  discs <- lapply(seq(1, 10, by = 0.5), function(x) disc(r = x))
  secondests <- secondorderprops(xiim, 
    gblargs = list(boxwidths = discs, estimators = c("GBLc", "GBLcc.pickaH", "GBLcc.mattfeldt", "GBLcc.pickaint", "GBLg.pickaint", "GBLg.mattfeldt")))
  expect_s3_class(secondests$GBL, "data.frame")
  expect_equal(nrow(secondests$GBL), length(discs))
  
  
  #test with discs and GBLemp - expect error
  discs <- lapply(seq(1, 10, by = 0.5), function(x) disc(r = x))
  expect_error(secondorderprops(xiim, 
    gblargs = list(boxwidths = discs, estimators = "GBLemp")),
    regexp = "non-numeric")
  
  #test with no centred covariance
  secondests <- secondorderprops(xiim, 
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

