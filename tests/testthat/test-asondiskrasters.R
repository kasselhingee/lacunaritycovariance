context("Format Conversion")

test_that("Conversion to and from an off RAM raster object works for simple lists", {
  ims <- lapply(heather, as.im, na.replace = 0)
  imsondisk <- imstoondiskras(ims)
  imsbackasim <- rasterstoims(imsondisk)
  
  #just tests first level structure of list
  expect_identical(summary(ims), summary(imsbackasim))
  
  expect_equivalent(lapply(imsondisk, inMemory), rep(FALSE, length(ims)))
  
  #expect the result to be identical up to unitname
  unitname(imsbackasim[[1]]) <- unitname(ims[[1]])
  unitname(imsbackasim[[2]]) <- unitname(ims[[2]])
  unitname(imsbackasim[[3]]) <- unitname(ims[[3]])
  expect_equal(ims, imsbackasim)
  #not quite identical() due to floating point accuracy of xrange and yrange.
})

test_that("Conversion to ondiskraster and back preserves structure of a complicated list", {
  ims <- lapply(heather, as.im, na.replace = 0)
  lsin <- c(list(im1 = ims), ims, miscval = "testging")
  rasondisk <- imstoondiskras(lsin, overwrite = TRUE)

  rasbacktoim <- rasterstoims(rasondisk)
  expect_equal(summary(lsin), summary(rasbacktoim))
})