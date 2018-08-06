context("Format Conversion")

test_that("getrastervaluesofpolygons gives an in-memory raster object", {
  suppressPackageStartupMessages(library("rgdal"))
  # rgdal is used here to read the ESRI shapefile into a SpatialPolygonsDataFrame
  regionfilepath <- system.file("extdata", package="stationaryracsinference")
  obspoly <- readOGR(regionfilepath, "aregionofinterest", verbose = FALSE)
  rsdatafilepath <- system.file("extdata/demorsraster.zip",
                              package="stationaryracsinference")
  rsfiles <- unzip(rsdatafilepath,exdir=tempdir())
  rasobj <- getrastervaluesofpolygons(obspoly, rsfiles[[2]])
  expect_false(fromDisk(rasobj[[1]]))
})
