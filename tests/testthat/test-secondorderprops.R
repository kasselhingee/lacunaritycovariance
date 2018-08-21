context("Estimation - second order properties")

test_that("secondorderprop() produces output of the expected classes given different arguments", {
  xi <- heather$coarse
  xiim <- as.im(xi, value = TRUE, na.replace = FALSE)

  #operates with all props going "everythingmode"
  expect_warning(secondests <- secondorderprops(xiim, 
    mvlargs = list(boxwidths = seq(1, 10, by = 0.1)),
    covarargs = list(estimators = "all"),
    paircorrargs = list(estimators = "all"),
    returnrotmean = FALSE),
    regexp = "harmonised")
  expect_s3_class(secondests$MVL, "fv")
  expect_s3_class(secondests$covariance, "imlist")
  expect_s3_class(secondests$paircorr, "imlist")
  expect_gt(length(secondests$MVL), 3)
  expect_gt(length(secondests$covariance), 3)
  expect_gt(length(secondests$paircorr), 3)

  #test with full defaults
  secondests <- secondorderprops(xiim)
  expect_is(secondests, "list")
  expect_length(secondests, 0)

  #test with rotmean and multiple estimators
  expect_warning(secondests <- secondorderprops(xiim, 
    mvlargs = list(boxwidths = seq(1, 10, by = 0.1)),
    covarargs = list(estimators = "all"),
    paircorrargs = list(estimators = "all"),
    returnrotmean = TRUE),
    regexp = "harmonised")
  expect_s3_class(secondests$MVL, "fv")
  expect_s3_class(secondests$covariance, "fv")
  expect_s3_class(secondests$paircorr, "fv")
  expect_gt(length(secondests$MVL), 3)
  expect_gt(length(secondests$covariance), 3)
  expect_gt(length(secondests$paircorr), 3)
  
  #test with rotmean and single estimators
  secondests <- secondorderprops(xiim, 
    mvlargs = list(boxwidths = seq(1, 10, by = 0.1), estimators = "MVLcc.pickaH"),
    covarargs = list(estimators = "pickaH"),
    paircorrargs = list(estimators = "pickaH"),
    returnrotmean = TRUE)
  expect_s3_class(secondests$MVL, "fv")
  expect_s3_class(secondests$covariance, "fv")
  expect_s3_class(secondests$paircorr, "fv")
  expect_length(secondests$MVL, 2)
  expect_length(secondests$covariance, 2)
  expect_length(secondests$paircorr, 2)
  
  #test without rotmean and single estimators
  secondests <- secondorderprops(xiim, 
    mvlargs = list(boxwidths = seq(1, 10, by = 0.1), estimators = "MVLcc.pickaH"),
    covarargs = list(estimators = "pickaH"),
    paircorrargs = list(estimators = "pickaH"),
    returnrotmean = FALSE)
  expect_s3_class(secondests$MVL, "fv")
  expect_s3_class(secondests$covariance, "imlist")
  expect_s3_class(secondests$paircorr, "imlist")
  expect_length(secondests$MVL, 2)
  expect_length(secondests$covariance, 1)
  expect_length(secondests$paircorr, 1)

  #test with only MVLgb and partially NULL for the other ones
  secondests <- secondorderprops(xiim, 
    mvlargs = list(boxwidths = seq(1, 10, by = 0.1), estimators = "MVLgb"))
  expect_s3_class(secondests$MVL, "fv")
  expect_length(secondests, 1)
  
  #test with discs and no MVLgb
  discs <- lapply(seq(1, 10, by = 0.5), function(x) disc(r = x))
  secondests <- secondorderprops(xiim, 
    mvlargs = list(boxwidths = discs, estimators = "MVLc"))
  expect_s3_class(secondests$MVL, "data.frame")
  expect_length(secondests, 1)
  
  #test with discs and multpile estimators,and no MVLgb
  discs <- lapply(seq(1, 10, by = 0.5), function(x) disc(r = x))
  secondests <- secondorderprops(xiim, 
    mvlargs = list(boxwidths = discs, estimators = c("MVLc", "MVLcc.pickaH", "MVLcc.mattfeldt", "MVLcc.pickaint", "MVLg.pickaint", "MVLg.mattfeldt")))
  expect_s3_class(secondests$MVL, "data.frame")
  expect_equal(nrow(secondests$MVL), length(discs))
  
  
  #test with discs and MVLgb - expect error
  discs <- lapply(seq(1, 10, by = 0.5), function(x) disc(r = x))
  expect_error(secondorderprops(xiim, 
    mvlargs = list(boxwidths = discs, estimators = "MVLgb")),
    regexp = "non-numeric")
})

