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

test_that("cppicka operates when passed im (logical and numeric) and fails correctly when im values not binary", {
  cpest <- cppicka(xiowin, obswin = w)
  expect_equal(max(abs(cpest - cppicka(xiim.l))), 0)
  expect_equal(max(abs(cpest - cppicka(xiim.n))), 0)
  expect_error(cppicka(xiim.n3), regexp = "has values other than")
})


test_that("racscovariance operates when passed im (logical and numeric) and fails correctly when im values not binary", {
  cvcest <- racscovariance(xiowin, obswin = w)
  expect_equal(max(abs(cvcest - racscovariance(xiim.l))), 0)
  expect_equal(max(abs(cvcest - racscovariance(xiim.n))), 0)
  expect_error(racscovariance(xiim.n3), regexp = "has values other than")
})

test_that("balancedracscovariance operates when passed im (logical and numeric) and fails correctly when im values not binary", {
  cvcest.o <- balancedracscovariances(xiowin, obswin = w)[["pickaH"]]
  cvcest.l <- balancedracscovariances(xiim.l)[["pickaH"]]
  cvcest.n <- balancedracscovariances(xiim.n)[["pickaH"]]
  expect_equal(max(abs(cvcest.o - cvcest.l)), 0)
  expect_equal(max(abs(cvcest.o - cvcest.n)), 0)
  expect_error(balancedracscovariances(xiim.n3), regexp = "has values other than")
})

test_that("pclns operates when passed im (logical and numeric) and fails correctly when im values not binary", {
  pclnest.o <- pclns(xiowin, obswin = w)[["pickaint"]]
  pclnest.l <- pclns(xiim.l)[["pickaint"]]
  pclnest.n <- pclns(xiim.n)[["pickaint"]]
  expect_equal(max(abs(pclnest.o - pclnest.l)), 0)
  expect_equal(max(abs(pclnest.o - pclnest.n)), 0)
  expect_error(pclns(xiim.n3), regexp = "has values other than")
})
