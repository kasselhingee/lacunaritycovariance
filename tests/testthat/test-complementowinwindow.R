context("Internals")

test_that("complement.owin.inwindow gives areas that are expected", {
  w <- disc(radius = 5, centre = c(5,10))
  xi <- heather$coarse
  xic <- complement.owin.inwindow(xi, w)
  expect_equal(area.owin(union.owin(intersect.owin(xi, w), xic)), area.owin(w), tolerance = area.owin(w)/1000)
  expect_equal(area.owin(intersect.owin(xi, xic)), 0)
  xip <- as.polygonal(xi)
  xipc <- complement.owin.inwindow(xip, w)
  expect_equal(area.owin(union.owin(intersect.owin(xip, w), xipc)), area.owin(w), tolerance = area.owin(w)/1000)
  expect_equal(area.owin(intersect.owin(xip, xipc)), 0, tolerance = area.owin(w)/1E6)
  xic <- complement.owin.inwindow(square(), square(r=0.5))
  expect_true(is.empty(xic))
  expect_equal(area.owin(xic), 0)
})
