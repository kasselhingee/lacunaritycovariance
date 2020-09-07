context("Estimation - GBL")

test_that("gblc() warns of unexpected inputs", {
  spatstat.options(npixel = 2^3)
  xiim_verytoy <- as.im(heather$coarse, na.replace = FALSE, eps = 2)
  covar <- plugincvc(xiim_verytoy)
  p <- sum(xiim_verytoy) * xiim_verytoy$xstep * xiim_verytoy$ystep / area(Frame(xiim_verytoy))
  sidelengths <- 2.2
  expect_error(gblc(sidelengths, covariance = covar, p = p, xiim = xiim_verytoy),
                 regexp = "Either covariance and p must be supplied or xiim supplied.")

  expect_error(gblc(sidelengths, covariance = covar, xiim = xiim_verytoy),
                 regexp = "Either covariance and p must be supplied or xiim supplied.")

  expect_error(gblc(sidelengths, p = p, xiim = xiim_verytoy),
                 regexp = "Either covariance and p must be supplied or xiim supplied.")
  reset.spatstat.options()
})

test_that("gblemp() warns of unexpected inputs", {
  sidel <- c(2.2)
  img <- as.im(heather$coarse,eps=c(2*heather$coarse$xstep, 4 * heather$coarse$xstep), na.replace = 0)
  expect_error(gblemp(sidel, img),
                 regexp = "image pixels must be square")

  expect_error(gblemp(sidel, 13),
                 regexp = "input xiim must be of class im")

})

test_that("GBLc estimates are historically consistent", {
  covar <- plugincvc(as.im(heather$coarse, na.replace = FALSE))
  p <- area(heather$coarse) / area(Frame(heather$coarse))
  sidelengths <- 2.2
  lac <- gblc(sidelengths, covar, p)
  expect_equal(lac$GBL, 1 + 0.05855459, tolerance = 0.02)
})

test_that("GBLemp estimates are historically consistent", {
  img <- as.im(heather$coarse, eps = heather$coarse$xstep, na.replace = 0)
  sidel <- c(2.2)
  lac <- gblemp(sidel, img)
  expect_equal(lac$GBL, 1 + 0.03253836)
})

test_that("GBLc estimates are consistent for input side lengths or owin squares", {
  spatstat.options(npixel = 2^3)
  xiim_verytoy <- as.im(heather$coarse, na.replace = FALSE, eps = 2)
  covar <- plugincvc(xiim_verytoy)
  p <- sum(xiim_verytoy) * xiim_verytoy$xstep * xiim_verytoy$ystep / area(Frame(xiim_verytoy))
  sidelengths <- 2.2
  lac <- gblc(sidelengths, covar, p)
  expect_equal(lac$GBL, gblc(list(square(2.2)), covar, p)$GBL)
  reset.spatstat.options()
})

test_that("GBLc estimates are operate on lists of owin objects", {
  spatstat.options(npixel = 2^3)
  discs <- lapply(seq(2, 10, by = 3), function(x) disc(r = x))
  lac <- gblc(discs, xiim = as.im(heather$coarse, na.replace = 0, eps = 2))
  expect_s3_class(lac, "data.frame")
  expect_equal(nrow(lac), length(discs))
  reset.spatstat.options()
})

test_that("Cubature integration with na.replace gives NA values when covariance argument is too large for covariance estimate", {
  skip_on_cran()
  spatstat.options(npixel = 2^3)
  xiim <- as.im(heather$coarse, na.replace = 0, eps = 2)
  suppressWarnings(ccov <- cencovariance(xi = xiim,
                        setcov_boundarythresh = 1E-20,
                        estimators = "pickaH"))
  p <- sum(xiim) / sum(is.finite(xiim$v))
  sidelengths <- c(8, 11)
  lac <- gblcc.inputcovar(sidelengths, ccov[[1]], p = p, integrationMethod = "cubature")
  expect_true(all(is.finite(lac$GBL) == c(TRUE, FALSE)))
  reset.spatstat.options()
})

test_that("harmonisesum integration with na.replace gives NA values when covariance argument is too large for covariance estimate", {
  spatstat.options(npixel = 2^3)
  sidelengths <- c(8, 11)
  xiim <- as.im(heather$coarse, na.replace = 0, eps = 2)
  suppressWarnings(ccov <- cencovariance(xi = xiim,
                        setcov_boundarythresh = 1E-20,
                        estimators = "pickaH"))
  p <- p <- sum(xiim) / sum(is.finite(xiim$v))
  lac <- gblcc.inputcovar(sidelengths, ccov[[1]], p = p, integrationMethod = "harmonisesum")
  expect_true(all(is.finite(lac$GBL) == c(TRUE, FALSE)))
  reset.spatstat.options()
})

test_that("gblcc estimates operate on lists of owin objects", {
  spatstat.options(npixel = 2^3)
  discs <- lapply(seq(2, 10, by = 3), function(x) disc(r = x))
  lac <- gblcc(discs, xiim = as.im(heather$coarse, na.replace = 0, eps = 2))
  expect_s3_class(lac, "data.frame")
  expect_equal(nrow(lac), length(discs))
  reset.spatstat.options()
})

test_that("gblg estimates operate on lists of owin objects", {
  spatstat.options(npixel = 2^3)
  discs <- lapply(seq(2, 10, by = 3), function(x) disc(r = x))
  lac <- gblg(discs, xiim = as.im(heather$coarse, na.replace = 0, eps = 2))
  expect_is(lac, "numeric")
  expect_length(lac, length(discs))
  reset.spatstat.options()
})

test_that("GBLc estimates are the same from estimated covariance or original image", {
  spatstat.options(npixel = 2^3)
  img <- as.im(heather$coarse, na.replace = 0, eps = 1)
  covar <- plugincvc(img)
  p <- sum(img) / sum(is.finite(img$v))
  sidelengths <- seq(1, 5, by = 2)
  gblc.covar <- gblc(sidelengths, covar, p)
  gblc.im <- gblc(sidelengths, xiim = img)
  expect_equal(gblc.covar, gblc.im)
  reset.spatstat.options()
})

test_that("integration when covar is constant gives squared area (i.e. gbl = 1)", {
  spatstat.options(npixel = 2^7)
  covar <- as.im(owin(c(-7, 7), c(-7, 7)), eps = 0.01)
  p <- 1
  sidelengths <- seq(1, 2.2, by = 0.5)
  lac <- gblc(sidelengths, covar, p)
  expect_equal(lac$GBL, rep(1, length(sidelengths)), tolerance = 0.01)
  
  expect_equal(gblc(lapply(c(0.5, 1, 3), disc), covar, p)$GBL, rep(1, 3), tolerance = 0.01)
  reset.spatstat.options()
})

test_that("GBLc and GBLemp produce similar results for large square observation windows", {
  skip_on_cran()
  #xiimg and covarest.frim is pregenerated in helper-calccovar
  sidelengths <- seq(xiimg$xstep * 3, 15, by = xiimg$xstep * 2) #odd pixel widths!
  lac.gblcc.picakHest <- gblc(sidelengths, covar = covarest.frim, p = xi.p)
  lac.gblempest <- gblemp(sidelengths, xiimg)

  expect_equal(lac.gblcc.picakHest$s, lac.gblempest$s)
  expect_equal(lac.gblcc.picakHest$GBL, lac.gblempest$GBL, tolerance = 5E-2)
})

test_that("gbl() fails nicely when GBLemp can't estimate anything", {
  spatstat.options(npixel = 2^2)
  xiim <- as.im(heather$coarse, value = TRUE, na.replace = FALSE, eps = 2)
  #fake lots of missing data
  xiim[shift.owin(reflect(heather$coarse), vec = c(10, 20))] <- NA
  expect_warning(gbl(xiim, seq(2, 10, by = 4), estimators = c("GBLcc.pickaH", "GBLemp")), regexp = "1 or fewer of the provided box widths")
  reset.spatstat.options()
})

test_that("gbl() harmonises estimates to produce meaningful fv object", {
  spatstat.options(npixel = 2^3)
  xiim_verytoy <- as.im(heather$coarse, value = TRUE, na.replace = FALSE, eps = 2)
  expect_warning(gblest <- gbl(xiim_verytoy, seq(2, 10, by = 3), estimators = c("GBLcc.pickaH", "GBLemp")), regexp = "harmon")
  expect_silent(lapply(gblest, plot.fv, limitsonly = TRUE))
  skip_on_cran()
  expect_silent(lapply(gblest, plot.fv, type = "n"))
  reset.spatstat.options()
})

test_that("gbl() operates nicely when only one estimator requested", {
  skip_on_cran() #less important test
  spatstat.options(npixel = 2^3)
  xiim <- as.im(heather$coarse, value = TRUE, na.replace = FALSE, eps = 1)
  expect_silent(gbl(xiim, seq(1, 10, by = 4), estimators = "GBLcc.pickaH"))
  reset.spatstat.options()
})

test_that("gbl() operates on owin style binary maps", {
  spatstat.options(npixel = 2^3)
  xiim <- as.im(heather$coarse, value = TRUE, na.replace = FALSE, eps = 1)
  xi <- as.mask(heather$coarse, eps = xiim$xstep)
  obswin <- setminus.owin(Frame(xi), square(5))
  xiim[square(5)] <- NA
  expect_warning(out <- gbl(xi,
                            seq(1, 10, by = 4),
                            estimators = c("GBLg.mattfeldt", "GBLg.pickaint", "GBLg.pickaH", "GBLcc.mattfeldt",
                                                                   "GBLcc.pickaint", "GBLc", "GBLemp"),
                            obswin = obswin))
  expect_warning(out_im <- gbl(xiim,
                               seq(1, 10, by = 4),
                               c("GBLg.mattfeldt", "GBLg.pickaint", "GBLg.pickaH", "GBLcc.mattfeldt",
                                 "GBLcc.pickaint", "GBLc", "GBLemp")))
  expect_equal(out, out_im)

  skip_on_cran() #less important tests mostly covered by above
  xiim <- as.im(heather$coarse, value = TRUE, na.replace = FALSE)
  xi <- heather$coarse
  obswin <- Frame(xi)
  expect_warning(out <- gbl(xi, 
                            seq(0.1, 10, by = 1),
                            estimators = c("GBLg.mattfeldt", "GBLg.pickaint", "GBLg.pickaH", "GBLcc.mattfeldt",
                                                                   "GBLcc.pickaint", "GBLc", "GBLemp"),
                            obswin = obswin))
  expect_warning(out_im <- gbl(xiim,
                               seq(0.1, 10, by = 1),
                            estimators = c("GBLg.mattfeldt", "GBLg.pickaint", "GBLg.pickaH", "GBLcc.mattfeldt",
                                                                   "GBLcc.pickaint", "GBLc", "GBLemp")
                               ))
  expect_equal(out, out_im)
  reset.spatstat.options()
})
