context("RACS simulation")

test_that("rbdr produces simulations with the correct coverage and covariance", {
  grain <- owin(xrange = c(-5, 5), yrange = c(-5, 5))
  win <- owin(xrange = c(0, 20), c(0, 20))
  if (identical(Sys.getenv("NOT_CRAN"), "true")){
    win <- owin(xrange = c(0, 200), c(0, 200))
  }
  lambda <- 4.2064E-3
  
  #calculate theoretical values of the model
  truecoveragefrac <- bdrcoverageprob(lambda, grain)
  xy <- as.mask(dilationAny(win, win), eps = c(0.2, 0.2))
   #eps of 1 (a 10 grain width) is too large! eps of 0.1
  thcovariance <- bdrcovar(lambda, grain, xy)
  ## use covariance to form approximate confidence interval
  setcovB <- setcov(win, eps = c(thcovariance$xstep, thcovariance$ystep))
  harmims <- harmonise.im(thcovar = thcovariance, setcovwin = setcovB)
  integrand <- eval.im((thcovariance - truecoveragefrac ^ 2) * setcovB, envir = list(thcovariance = harmims$thcovar, setcovB = harmims$setcovwin))
  exactvariance <- (1 / (area.owin(win))^2) * sum(integrand) * integrand$xstep * integrand$ystep
  ##assume that window large enough that distribution of estimator is Gaussian
  ##So half-width of 0.05 conf interval is
  confint_halfwidth <- -qnorm(0.025, mean = 0, sd = sqrt(exactvariance))

  xi <- rbdr(lambda, grain, win)
  phat <- coverageprob(xi, win)
  # plot(xi, col = "black")
  # plot(win, add = TRUE)
  expect_equal(phat, truecoveragefrac, tolerance = confint_halfwidth)
  
  #now on to checking covariance
  cvchat <- racscovariance(as.mask(xi, eps = c(0.2, 0.2)), obswin = win, estimators = list("pickaint"), drop = TRUE)
  truecvc.iso <- rotmean(thcovariance[disc(radius = 50), drop = FALSE], padzero = FALSE)
  cvchat.iso <- rotmean(cvchat[disc(radius = 50), drop = FALSE], padzero = FALSE)
  cvciso <- collapse.fv(harmonise.fv(truecvc.iso, cvchat.iso), different = "f")
  expect_lt(with.fv(cvciso, max(x1 - x2)), 10 * confint_halfwidth) #the 10 is just a guess based on 4th order properties being much harder to estimate
})
