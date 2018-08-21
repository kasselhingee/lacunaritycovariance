context("RACS simulation")

test_that("rbdd produces simulations with the correct area fraction", {
  #Boolean model with discs of radius 10.
  #The intensity has been chosen such that the true coverage fraction is very close to 0.5.
  discr <- 10
  w <- owin(xrange = c(0, 100), c(0, 100))
  lambda <- 2.2064E-3
  xi <- rbdd(lambda, discr, w, seed=6549)
  expect_equal(coveragefrac(xi, w), 0.4644903, tolerance = 1E-7)

  #theoretical coveragefrac
  truecoveragefrac <- bddcoverageprob(lambda, discr)

  ## use covariance to form approximate confidence interval
  thcovariance <- bddcovar(xrange = c(-100, 100), yrange = c(-100, 100), eps = c(1, 1), lambda, discr)
  setcovB <- setcov(w, eps = c(thcovariance$xstep, thcovariance$ystep))
  harmims <- harmonise.im(thcovar = thcovariance, setcovwin = setcovB)
  integrand <- eval.im((thcovariance - truecoveragefrac ^ 2) * setcovB, envir = list(thcovariance = harmims$thcovar, setcovB = harmims$setcovwin))
  exactvariance <- (1 / (area.owin(w))^2) * sum(integrand) * integrand$xstep * integrand$ystep
  ##assume that window large enough that distribution of estimator is Gaussian
  ##So half-width of 0.05 conf interval is
  confint_halfwidth <- -qnorm(0.025, mean = 0, sd = sqrt(exactvariance))

  #test on random simulations (seed unfixed)
  xi <- rbdd(lambda, discr, w)
  expect_equal(coveragefrac(xi, w), truecoveragefrac, tolerance = confint_halfwidth)
})

test_that("Computations of germ intensity, disc radius and coverage probability are correct",
 {
  expect_equal(bddcoverageprob(0.005, 5), 0.3247681)
  expect_equal(bddlambda(bddcoverageprob(0.005, 5), 5), 0.005)
  expect_equal(bdddiscr(bddcoverageprob(0.005, 5), 0.005), 5)
})
