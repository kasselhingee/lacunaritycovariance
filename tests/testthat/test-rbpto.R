context("RACS Simulation")

test_that("rbpto generates simulations that match covariance for big window", {
  win <- owin()
  grain <- disc(r = 0.01)
  lambda <- 500
  xm <- 0.01
  alpha <- 2
  lengthscales <- seq(1, 5, by = 0.1)
  cp <- bpto.coverageprob(lambda, grain, xm, alpha, lengthscales = lengthscales)
  xis <- replicate(3, rbpto(lambda, grain, win, xm, alpha,
	      lengthscales = lengthscales), simplify = FALSE)

  #compare coverage probability to phat
  expect_equal(mean(vapply(xis, coverageprob, obswin = win, FUN.VALUE = 0.0)),
  cp, tol = 5E-2)

  #compare covariance to estimated covariance
  xy <- as.mask(win, eps = 0.001)
  cvc.th <- bpto.covar(lambda, grain, xm, alpha, lengthscales = lengthscales, xy)
  xiims <- lapply(xis, as.im.owin, W = win, eps = 0.001, value = TRUE, na.replace = FALSE)
  cvc.ests <- lapply(xiims, function(x) racscovariance(x, estimators = "pickaH", drop = TRUE))
  cvc.est <- summary.imlist(cvc.ests)$mean

  #compare to cvc.th
  expect_warning(cvc.diff <- eval.im(cvc.th - cvc.est, harmonize = TRUE), regexp = "images .* were not compatible")
  #plot(solist(cvc.th, cvc.est, cvc.diff[Frame(cvc.th)]), axes = TRUE)
  # the differences make very beautiful patterns
  expect_lt(max(abs(cvc.diff)), cp * 1E-1)
  expect_lt(mean(abs(cvc.diff)), cp * 5E-2)
})

test_that("bpto.coverageprob is consistent with bpto.germintensity", {
  grain <- disc(r = 0.01)
  lambda <- 20
  xm <- 0.01
  alpha <- 2
  lengthscales <- seq(1, 5, by = 0.1)
  expect_equal(lambda, bpto.germintensity( bpto.coverageprob(lambda, grain, xm, alpha, lengthscales = lengthscales),
      grain, xm, alpha, lengthscales = lengthscales) )
})
