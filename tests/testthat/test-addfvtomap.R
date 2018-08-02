context("Plotting Functions")

test_that("addfvtomap doesn't error", {
  plot(heather$coarse, axes = TRUE)
  fvobj <- Hest(heather$coarse)
  expect_silent(addfvtomap(fvobj, mapxlim = c(0, 10), mapylim = c(5, 10)))
})