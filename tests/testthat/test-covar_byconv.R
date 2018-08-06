context("Covariance Estimation")

test_that("balanced covar estimation using convolutions as the base object match balancedracscovariance", {
  out1 <- balancedracscovariances(heather$coarse, obswin = Frame(heather$coarse))
  out2 <- byconv_cvchats(heather$coarse, obswin = Frame(heather$coarse))
  expect_equal(names(out1), names(out2))
  suppressWarnings(diffs <- mapply(function(x, y) eval.im(x - y), out1, out2, SIMPLIFY = FALSE))
  meandiff <- vapply(diffs, function(x) mean(x), FUN.VALUE = 0.0)
  expect_lt(max(meandiff), phat^2/100)
  meansqdiff <- vapply(diffs, function(x) mean(x^2), FUN.VALUE = 0.0)
  expect_lt(max(meansqdiff), (phat*10^(-2))^2)
  mdiffs <- vapply(diffs, function(x) max(abs(x)), FUN.VALUE = 0.0)
  expect_lt(max(mdiffs), phat^2/10)

  #Note that above results are slightly different to due to byconv <- cvchats accepting an xy argument of the as.mask of windows for setcov.
  #the following should skip this difference

  xi <- heather$coarse
  obswin <- Frame(xi)
  xixi <- setcov(xi, xy = xi)
  winwin <- setcov(obswin, eps = c(xixi$xstep, xixi$ystep))
  xiwin <- setcov(xi, obswin, xy = xi)
  phat <- coverageprob(xi, obswin = Frame(xi))

  out3 <- cvchats_convolves(xixi, winwin, xiwin = xiwin, phat = phat, modifications = "all")
  diffs <- mapply(function(x, y) eval.im(x - y), out1, out3, SIMPLIFY = FALSE)
  meandiff <- vapply(diffs, function(x) mean(x), FUN.VALUE = 0.0)
  expect_lt(max(meandiff), phat^2/100)
  meansqdiff <- vapply(diffs, function(x) mean(x^2), FUN.VALUE = 0.0)
  expect_lt(max(meansqdiff), (phat*10^(-2))^2)
  mdiffs <- vapply(diffs, function(x) max(abs(x)), FUN.VALUE = 0.0)
  expect_lt(max(mdiffs), phat^2/10)
})
