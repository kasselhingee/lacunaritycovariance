context("Internals")

spatstat.options(npixel = 2^3)

#creating a test that involves NA values in an image
xifull <- as.mask(heather$coarse, eps = 5 * heather$coarse$xstep)
w <- setminus.owin(Frame(xifull), square(5)) 
xiowin <- intersect.owin(xifull, w)
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


test_that("plugincvc operates when passed im (logical and numeric) and fails correctly when im values not binary", {
  cvcest <- plugincvc(xiowin, obswin = w)
  expect_equal(max(abs(cvcest - plugincvc(xiim.l))), 0)
  expect_equal(max(abs(cvcest - plugincvc(xiim.n))), 0)
  expect_error(plugincvc(xiim.n3), regexp = "has values other than")
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

test_that("gblemp operates as expected when pass im object", {
  boxwidths <- seq(1, 5, by = 1)
  gblemp.l <- gblemp(boxwidths, xiim.l)
  gblemp.n <- gblemp(boxwidths, xiim.n)
  expect_equal(max(abs(gblemp.l - gblemp.n)), 0)
  expect_error(gblemp(boxwidths, xiim.n3), regexp = "has values other than")
})

test_that("gbl() operates well when passed im or owin boxes", {
  boxwidths <- seq(1, 5, by = 1)
  expect_warning(gblest.l <- gbl(xiim.l, boxwidths,
                                 c("GBLg.mattfeldt", "GBLg.pickaint", "GBLg.pickaH", "GBLcc.mattfeldt",
                                   "GBLcc.pickaint", "GBLc", "GBLemp")),
                 regexp = "harmon")
  expect_warning(gblest.n <- gbl(xiim.n, boxwidths,
                                 c("GBLg.mattfeldt", "GBLg.pickaint", "GBLg.pickaH", "GBLcc.mattfeldt",
                                   "GBLcc.pickaint", "GBLc", "GBLemp")),
                 regexp = "harmon")
  expect_equal(gblest.l, gblest.n)
  expect_error(gbl(xiim.n3, boxwidths), regexp = "isbinarymap")
})

reset.spatstat.options()
