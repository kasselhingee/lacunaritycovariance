context("Estimation - Contagion")

test_that("contagdiscstate is returning the correct format", {
  xi <- heather$coarse
  obswindow <- Frame(heather$coarse)
  p <- coverageprob(xi, Frame(xi))
  xiH <- Hest(xi, W = obswindow) #Sph. Contact Distrution Estimate
  xicH <- Hest(complement.owin(xi), W = obswindow) #Conditional Core Prob. Estimate
  contagion <- contagdiscstate(xiH, xicH, p, normalise = TRUE)

  expect_is(contagion, "fv")
  expect_equal(names(contagion), as.character(list("r", "contag")))

  coldescriptions <- attr(contagion, "desc")
  expect_equal(coldescriptions, as.character(list("radius", "normalised SCD contagion estimate")))

  contagion <- contagdiscstate(xiH, xicH, p, normalise = FALSE)
  coldescriptions <- attr(contagion, "desc")
  expect_equal(coldescriptions, as.character(list("radius", "unnormalised SCD contagion estimate")))
})
