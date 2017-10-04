context("MVL estimation")

test_that("mvlc() warns of unexpected inputs", {
  covar <- racscovariance(heather$coarse)
  p <- area(heather$coarse) / area(Frame(heather$coarse))
  sidelengths <- 2.2
  img <- as.im(heather$coarse,eps=heather$coarse$xstep, na.replace=0)
  expect_error(lac(sidelengths, covariance = covar, p = p, xiim = img),
                 regexp = "Either covariance and p must be supplied or xiim supplied.")
  
  expect_error(lac(sidelengths, covariance = covar, xiim = img),
                 regexp = "Either covariance and p must be supplied or xiim supplied.")
  
  expect_error(lac(sidelengths, p = p, xiim = img),
                 regexp = "Either covariance and p must be supplied or xiim supplied.")
})

test_that("MVLc estimates are historically consistent", {
  covar <- racscovariance(heather$coarse)
  p <- area(heather$coarse) / area(Frame(heather$coarse))
  sidelengths <- 2.2
  lac <- lac(sidelengths, covar, p)
  expect_equal(lac$MVL, 0.05855459, tolerance = 1E-6)
})

test_that("MVLgb estimates are historically consistent", {
  img <- as.im(heather$coarse,eps=heather$coarse$xstep, na.replace=0)
  sidel <- c(2.2)
  lac.wraw <- lacgb(img,sidel,inclraw = TRUE)
  expect_equal(lac.wraw$MVL, 0.03253836)
  expect_equal(lac.wraw$raw, -0.05775767)
})

test_that("MVLc estimates are consitent for input side lengths or owin squares", {
  covar <- racscovariance(heather$coarse)
  p <- area(heather$coarse) / area(Frame(heather$coarse))
  sidelengths <- 2.2
  lac <- lac(sidelengths, covar, p)
  expect_equal(lac$MVL, mvlc(list(square(2.2)), covar, p))
})

test_that("integration when covar is constant gives squared area", {
  covar <- as.im(owin(c(-6, 6), c(-6, 6)), eps = 0.01)
  p <- 1
  sidelengths <- seq(1, 2.2, by = 0.1)
  lac <- lac(sidelengths, covar, p)
  expect_equal(lac$MVL, rep(0, length(sidelengths)), tolerance = 0.01)
  
  expect_equal(lac(lapply(c(0.5, 1, 2, 3), disc), covar, p), rep(0, 4), tolerance = 0.01)
})

test_that("MVLc and MVLgb produce similar results for large square observation windows", {
  lambda <- 4 * 2.2064E-3
  discr <- 5
  w <- owin(xrange = c(0, 100) * 3, yrange = c(0, 100) * 3)
  xi <- rbdd(lambda, discr,w)
  xiimg <- as.im(xi, W = w, eps = c(0.1, 0.1), na.replace = 0)
  
  sidelengths = seq(0.3, 15, by = 0.2) #odd pixel widths!
  lac.mvlcest <- mvlc(sidelengths, xiim = xiimg)
  lac.mvlgbest <- mvlgb(xiimg, sidelengths)
  
  expect_equal(lac.mvlcest$s, lac.mvlgbest$s)
  expect_equal(lac.mvlcest$MVL, lac.mvlgbest$MVL, tolerance = 1E-2)
})
