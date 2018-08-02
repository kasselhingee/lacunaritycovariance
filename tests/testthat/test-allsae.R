context("Estimating Cover Area")

test_that("allsae() produces plausible output", {
  xi <- heather$coarse
  obswin <- Frame(heather$coarse)
  corrrad <- 2
  corrstepheight <- 0.95
  n11 <- 90
  n21 <- 10
  n22 <- 95
  n12 <- 5
  expect_warning(estimates <- allsae(xi, obswin, corrrad, corrstepheight, erosionrad = 5, n11, n21, n12, n22),
                 regexp = "owin converted to image using default pixel amounts")
  expect_equal(
    estimates[rownames(estimates) != "ErodeDilate", "areahat"],
    rep(estimates[1, 1], sum(rownames(estimates) != "ErodeDilate")))
  expect_equivalent(dim(estimates), c(8, 2))
})