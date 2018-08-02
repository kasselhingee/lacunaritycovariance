context("Plotting Functions")

test_that("addfvtomap doesn't error", {
  fvobj <- Hest(heather$coarse)
  plot(heather$coarse, axes = TRUE)
  expect_silent(addfvtomap(fvobj, mapxlim = c(0, 10), mapylim = c(5, 10)))
  if(require(vdiffr)){
    expect_doppelganger("map with fv on top", {
      plot(heather$coarse)
      addfvtomap(fvobj, mapxlim = c(0, 10), mapylim = c(5, 10))
    })
  }
})