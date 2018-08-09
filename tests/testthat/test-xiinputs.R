context("Reading Ims")

#creating a test that involves NA values in an image
w <- setminus.owin(Frame(heather$coarse), square(5)) 
xiowin <- intersect.owin(heather$coarse, w)
phat <- coverageprob(xiowin, obswin = w)
  
#creating logically valued image to pass
xiim.l <- as.im(xiowin, W = w, value = TRUE, na.replace = FALSE)
xiim.l[ setminus.owin(Frame(xiim.l), w) ] <- NA

#numerically valued image of 0, 1 and NA.
xiim.n <- xiim.l * 1

#numerically valued image of 0, 3 and NA.
xiim.n3 <- xiim.l * 3

test_that("isbinarymap operates correctly", {
  expect_true(isbinarymap(xiim.n))
  expect_true(isbinarymap(xiim.l))
  expect_false(isbinarymap(xiim.n3))
  expect_error(isbinarymap(xiim.n3, requiretrue = TRUE))
})

test_that("cppicka operates when passed im (logical and numeric) and fails correctly when im values not binary", {
  cpest <- cppicka(xiowin, obswin = w)
  expect_equal(max(abs(cpest - cppicka(xiim.l))), 0)
  expect_equal(max(abs(cpest - cppicka(xiim.n))), 0)
  expect_error(cppicka(xiim.n3), regexp = "has values other than")
})


test_that("tradcovarest operates when passed im (logical and numeric) and fails correctly when im values not binary", {
  cvcest <- tradcovarest(xiowin, obswin = w)
  expect_equal(max(abs(cvcest - tradcovarest(xiim.l))), 0)
  expect_equal(max(abs(cvcest - tradcovarest(xiim.n))), 0)
  expect_error(tradcovarest(xiim.n3), regexp = "has values other than")
})

test_that("balancedracscovariance operates when passed im (logical and numeric) and fails correctly when im values not binary", {
  cvcest.o <- racscovariance(xiowin, obswin = w)[["pickaH"]]
  cvcest.l <- racscovariance(xiim.l)[["pickaH"]]
  cvcest.n <- racscovariance(xiim.n)[["pickaH"]]
  expect_equal(max(abs(cvcest.o - cvcest.l)), 0)
  expect_equal(max(abs(cvcest.o - cvcest.n)), 0)
  expect_error(racscovariance(xiim.n3), regexp = "has values other than")
})

test_that("paircorr operates when passed im (logical and numeric) and fails correctly when im values not binary", {
  pclnest.o <- paircorr(xiowin, obswin = w)[["pickaint"]]
  pclnest.l <- paircorr(xiim.l)[["pickaint"]]
  pclnest.n <- paircorr(xiim.n)[["pickaint"]]
  expect_equal(max(abs(pclnest.o - pclnest.l)), 0)
  expect_equal(max(abs(pclnest.o - pclnest.n)), 0)
  expect_error(paircorr(xiim.n3), regexp = "has values other than")
})

test_that("mvlgb operates as expected when pass im object", {
  boxwidths <- seq(1, 5, by = 1)
  mvlgb.l <- mvlgb(boxwidths, xiim.l)
  mvlgb.n <- mvlgb(boxwidths, xiim.n)
  expect_equal(max(abs(mvlgb.l - mvlgb.n)), 0)
  expect_error(mvlgb(boxwidths, xiim.n3), regexp = "has values other than")
})

test_that("mvl() operates well when passed im or owin objects", {
  boxwidths <- seq(1, 5, by = 1)
  expect_warning(mvlest.l <- mvl(xiim.l, boxwidths), regexp = "harmon")
  expect_warning(mvlest.n <- mvl(xiim.n, boxwidths), regexp = "harmon")
  expect_equal(mvlest.l, mvlest.n)
  expect_error(mvl(xiim.n3, boxwidths), regexp = "has values other than")
})
