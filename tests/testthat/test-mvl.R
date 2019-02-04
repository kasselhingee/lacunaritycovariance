context("Estimation - GBL")

test_that("gblc() warns of unexpected inputs", {
  covar <- tradcovarest(heather$coarse)
  p <- area(heather$coarse) / area(Frame(heather$coarse))
  sidelengths <- 2.2
  img <- as.im(heather$coarse, eps = heather$coarse$xstep, na.replace = 0)
  expect_error(gblc(sidelengths, covariance = covar, p = p, xiim = img),
                 regexp = "Either covariance and p must be supplied or xiim supplied.")

  expect_error(gblc(sidelengths, covariance = covar, xiim = img),
                 regexp = "Either covariance and p must be supplied or xiim supplied.")

  expect_error(gblc(sidelengths, p = p, xiim = img),
                 regexp = "Either covariance and p must be supplied or xiim supplied.")
})

test_that("gblemp() warns of unexpected inputs", {
  sidel <- c(2.2)
  img <- as.im(heather$coarse,eps=c(heather$coarse$xstep, 2 * heather$coarse$xstep), na.replace = 0)
  expect_error(gblemp(sidel, img),
                 regexp = "image pixels must be square")

  expect_error(gblemp(sidel, 13),
                 regexp = "input xiim must be of class im")

})

test_that("GBLc estimates are historically consistent", {
  covar <- tradcovarest(heather$coarse)
  p <- area(heather$coarse) / area(Frame(heather$coarse))
  sidelengths <- 2.2
  lac <- gblc(sidelengths, covar, p)
  expect_equal(lac$GBL, 1 + 0.05855459, tolerance = 1E-6)
})

test_that("GBLgb estimates are historically consistent", {
  img <- as.im(heather$coarse, eps = heather$coarse$xstep, na.replace = 0)
  sidel <- c(2.2)
  lac <- gblemp(sidel, img)
  expect_equal(lac$GBL, 1 + 0.03253836)
})

test_that("GBLc estimates are consistent for input side lengths or owin squares", {
  covar <- tradcovarest(heather$coarse)
  p <- area(heather$coarse) / area(Frame(heather$coarse))
  sidelengths <- 2.2
  lac <- gblc(sidelengths, covar, p)
  expect_equal(lac$GBL, gblc(list(square(2.2)), covar, p)$GBL)
})

test_that("GBLc estimates are operate on lists of owin objects", {
  discs <- lapply(seq(1, 10, by = 0.5), function(x) disc(r = x))
  lac <- gblc(discs, xiim = as.im(heather$coarse, na.replace = 0))
  expect_s3_class(lac, "data.frame")
  expect_equal(nrow(lac), length(discs))
})

test_that("gblcc estimates operate on lists of owin objects", {
  discs <- lapply(seq(1, 10, by = 0.5), function(x) disc(r = x))
  lac <- gblcc(discs, xiim = as.im(heather$coarse, na.replace = 0))
  expect_s3_class(lac, "data.frame")
  expect_equal(nrow(lac), length(discs))
})

test_that("gblg estimates operate on lists of owin objects", {
  discs <- lapply(seq(1, 10, by = 0.5), function(x) disc(r = x))
  lac <- gblg(discs, xiim = as.im(heather$coarse, na.replace = 0))
  expect_is(lac, "numeric")
  expect_length(lac, length(discs))
})

test_that("GBLc estimates are the same from estimated covariance or original image", {
  img <- as.im(heather$coarse, eps = heather$coarse$xstep, na.replace = 0)
  covar <- tradcovarest(img)
  p <- sum(img) / sum(is.finite(img$v))
  sidelengths <- seq(1, 5, by = heather$coarse$xstep*2)
  gblc.covar <- gblc(sidelengths, covar, p)
  gblc.im <- gblc(sidelengths, xiim = img)
  expect_equal(gblc.covar, gblc.im)
})

test_that("integration when covar is constant gives squared area (i.e. gbl = 1)", {
  covar <- as.im(owin(c(-6, 6), c(-6, 6)), eps = 0.01)
  p <- 1
  sidelengths <- seq(1, 2.2, by = 0.1)
  lac <- gblc(sidelengths, covar, p)
  expect_equal(lac$GBL, rep(1, length(sidelengths)), tolerance = 0.01)
  
  expect_equal(gblc(lapply(c(0.5, 1, 2, 3), disc), covar, p)$GBL, rep(1, 4), tolerance = 0.01)
})

test_that("GBLc and GBLgb produce similar results for large square observation windows", {
  #xiimg and covarest.frim is pregenerated in helper-calccovar
  sidelengths <- seq(xiimg$xstep * 3, 15, by = xiimg$xstep * 2) #odd pixel widths!
  lac.gblcc.picakHest <- gblc(sidelengths, covar = covarest.frim, p = xi.p)
  lac.gblempest <- gblemp(sidelengths, xiimg)

  expect_equal(lac.gblcc.picakHest$s, lac.gblempest$s)
  expect_equal(lac.gblcc.picakHest$GBL, lac.gblempest$GBL, tolerance = 5E-2)
})

test_that("gbl() fails nicely when GBLgb can't estimate anything", {
  xiim <- as.im(heather$coarse, value = TRUE, na.replace = FALSE)
  #fake lots of missing data
  xiim[shift.owin(reflect(heather$coarse), vec = c(10, 20))] <- NA
  expect_warning(gbl(xiim, seq(1, 10, by = 1)), regexp = "1 or fewer of the provided box widths")
})

test_that("gbl() harmonises estimates to produce meaningful fv object", {
  xiim <- as.im(heather$coarse, value = TRUE, na.replace = FALSE)
  expect_warning(gblest <- gbl(xiim, seq(0.2, 10, by = 1)), regexp = "harmon")
  expect_silent(lapply(gblest, plot.fv, type = "n"))
})

test_that("gbl() operates nicely when only one estimator requested", {
  xiim <- as.im(heather$coarse, value = TRUE, na.replace = FALSE)
  expect_silent(gbl(xiim, seq(0.1, 10, by = 1), estimators = "GBLc"))
})

test_that("gbl() operates on owin style binary maps", {
  xiim <- as.im(heather$coarse, value = TRUE, na.replace = FALSE)
  xi <- heather$coarse
  obswin <- setminus.owin(Frame(heather$coarse), square(5))
  xiim[square(5)] <- NA
  expect_warning(out <- gbl(xi, seq(0.1, 10, by = 1), obswin = obswin))
  expect_warning(out_im <- gbl(xiim, seq(0.1, 10, by = 1)))
  expect_equal(out, out_im)
  
  xiim <- as.im(heather$coarse, value = TRUE, na.replace = FALSE)
  xi <- heather$coarse
  obswin <- Frame(xi)
  expect_warning(out <- gbl(xi, seq(0.1, 10, by = 1), obswin = obswin))
  expect_warning(out_im <- gbl(xiim, seq(0.1, 10, by = 1)))
  expect_equal(out, out_im)
})
