context("Estimation - MVL")

test_that("mvlc() warns of unexpected inputs", {
  covar <- tradcovarest(heather$coarse)
  p <- area(heather$coarse) / area(Frame(heather$coarse))
  sidelengths <- 2.2
  img <- as.im(heather$coarse, eps = heather$coarse$xstep, na.replace = 0)
  expect_error(mvlc(sidelengths, covariance = covar, p = p, xiim = img),
                 regexp = "Either covariance and p must be supplied or xiim supplied.")

  expect_error(mvlc(sidelengths, covariance = covar, xiim = img),
                 regexp = "Either covariance and p must be supplied or xiim supplied.")

  expect_error(mvlc(sidelengths, p = p, xiim = img),
                 regexp = "Either covariance and p must be supplied or xiim supplied.")
})

test_that("mvlgb() warns of unexpected inputs", {
  sidel <- c(2.2)
  img <- as.im(heather$coarse,eps=c(heather$coarse$xstep, 2 * heather$coarse$xstep), na.replace = 0)
  expect_error(mvlgb(sidel, img),
                 regexp = "image pixels must be square")

  expect_error(mvlgb(sidel, 13),
                 regexp = "input xiim must be of class im")

})

test_that("MVLc estimates are historically consistent", {
  covar <- tradcovarest(heather$coarse)
  p <- area(heather$coarse) / area(Frame(heather$coarse))
  sidelengths <- 2.2
  lac <- mvlc(sidelengths, covar, p)
  expect_equal(lac$MVL, 0.05855459, tolerance = 1E-6)
})

test_that("MVLgb estimates are historically consistent", {
  img <- as.im(heather$coarse, eps = heather$coarse$xstep, na.replace = 0)
  sidel <- c(2.2)
  lac <- mvlgb(sidel, img)
  expect_equal(lac$MVL, 0.03253836)
})

test_that("MVLc estimates are consistent for input side lengths or owin squares", {
  covar <- tradcovarest(heather$coarse)
  p <- area(heather$coarse) / area(Frame(heather$coarse))
  sidelengths <- 2.2
  lac <- mvlc(sidelengths, covar, p)
  expect_equal(lac$MVL, mvlc(list(square(2.2)), covar, p)$MVL)
})

test_that("MVLc estimates are operate on lists of owin objects", {
  discs <- lapply(seq(1, 10, by = 0.5), function(x) disc(r = x))
  lac <- mvlc(discs, xiim = as.im(heather$coarse, na.replace = 0))
  expect_s3_class(lac, "data.frame")
  expect_equal(nrow(lac), length(discs))
})

test_that("mvlcc estimates operate on lists of owin objects", {
  discs <- lapply(seq(1, 10, by = 0.5), function(x) disc(r = x))
  lac <- mvlcc(discs, xiim = as.im(heather$coarse, na.replace = 0))
  expect_s3_class(lac, "data.frame")
  expect_equal(nrow(lac), length(discs))
})

test_that("mvlg estimates operate on lists of owin objects", {
  discs <- lapply(seq(1, 10, by = 0.5), function(x) disc(r = x))
  lac <- mvlg(discs, xiim = as.im(heather$coarse, na.replace = 0))
  expect_is(lac, "numeric")
  expect_length(lac, length(discs))
})

test_that("MVLc estimates are the same from estimated covariance or original image", {
  img <- as.im(heather$coarse, eps = heather$coarse$xstep, na.replace = 0)
  covar <- tradcovarest(img)
  p <- sum(img) / sum(is.finite(img$v))
  sidelengths <- seq(1, 5, by = heather$coarse$xstep*2)
  mvlc.covar <- mvlc(sidelengths, covar, p)
  mvlc.im <- mvlc(sidelengths, xiim = img)
  expect_equal(mvlc.covar, mvlc.im)
})

test_that("integration when covar is constant gives squared area", {
  covar <- as.im(owin(c(-6, 6), c(-6, 6)), eps = 0.01)
  p <- 1
  sidelengths <- seq(1, 2.2, by = 0.1)
  lac <- mvlc(sidelengths, covar, p)
  expect_equal(lac$MVL, rep(0, length(sidelengths)), tolerance = 0.01)
  
  expect_equal(mvlc(lapply(c(0.5, 1, 2, 3), disc), covar, p)$MVL, rep(0, 4), tolerance = 0.01)
})

test_that("MVLc and MVLgb produce similar results for large square observation windows", {
  #xiimg and covarest.frim is pregenerated in helper-calccovar
  sidelengths <- seq(xiimg$xstep * 3, 15, by = xiimg$xstep * 2) #odd pixel widths!
  lac.mvlcc.picakHest <- mvlc(sidelengths, covar = covarest.frim, p = xi.p)
  lac.mvlgbest <- mvlgb(sidelengths, xiimg)

  expect_equal(lac.mvlcc.picakHest$s, lac.mvlgbest$s)
  expect_equal(lac.mvlcc.picakHest$MVL, lac.mvlgbest$MVL, tolerance = 5E-2)
})

test_that("mvl() fails nicely when MVLgb can't estimate anything", {
  xiim <- as.im(heather$coarse, value = TRUE, na.replace = FALSE)
  #fake lots of missing data
  xiim[shift.owin(reflect(heather$coarse), vec = c(10, 20))] <- NA
  expect_warning(mvl(xiim, seq(1, 10, by = 1)), regexp = "1 or fewer of the provided box widths")
})

test_that("mvl() harmonises estimates to produce meaningful fv object", {
  xiim <- as.im(heather$coarse, value = TRUE, na.replace = FALSE)
  expect_warning(mvlest <- mvl(xiim, seq(0.2, 10, by = 1)), regexp = "harmon")
  expect_silent(lapply(mvlest, plot.fv, type = "n"))
})

test_that("mvl() operates nicely when only one estimator requested", {
  xiim <- as.im(heather$coarse, value = TRUE, na.replace = FALSE)
  expect_silent(mvl(xiim, seq(0.1, 10, by = 1), estimators = "MVLc"))
})

test_that("mvl() operates on owin style binary maps", {
  xiim <- as.im(heather$coarse, value = TRUE, na.replace = FALSE)
  xi <- heather$coarse
  obswin <- setminus.owin(Frame(heather$coarse), square(5))
  xiim[square(5)] <- NA
  expect_warning(out <- mvl(xi, seq(0.1, 10, by = 1), obswin = obswin))
  expect_warning(out_im <- mvl(xiim, seq(0.1, 10, by = 1)))
  expect_equal(out, out_im)
  
  xiim <- as.im(heather$coarse, value = TRUE, na.replace = FALSE)
  xi <- heather$coarse
  obswin <- Frame(xi)
  expect_warning(out <- mvl(xi, seq(0.1, 10, by = 1), obswin = obswin))
  expect_warning(out_im <- mvl(xiim, seq(0.1, 10, by = 1)))
  expect_equal(out, out_im)
})
